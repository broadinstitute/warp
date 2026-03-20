"""
Workflow runner for test harness.

Orchestrates the complete test execution workflow including Docker setup,
GCS data preparation, WDL patching, and miniwdl execution.
"""

import subprocess
import os
from pathlib import Path
from typing import Optional, Tuple
from .config import TestConfig
from .docker_manager import DockerManager
from .gcs_manager import GCSManager
from .wdl_patcher import WDLPatcher, create_simple_docker_replacements


class WorkflowRunner:
    """Orchestrates workflow test execution."""
    
    def __init__(self, config: TestConfig):
        """Initialize runner.
        
        Args:
            config: TestConfig with all necessary settings
        """
        self.config = config
        self.docker_manager = DockerManager(config.docker)
        self.gcs_manager = GCSManager(config.gcs, config.repo_root)
        self.wdl_patcher = WDLPatcher(config.repo_root)
    
    def cleanup_old_runs(self) -> bool:
        """Clean up old miniwdl run directories.
        
        Removes directories matching pattern: YYYYMMDD_HHMMSS_<workflow_name>
        
        Returns:
            True if successful
        """
        if self.config.skip_cleanup:
            print("Skipping cleanup (SKIP_CLEANUP=true)")
            return True
        
        print("Cleaning up old miniwdl runs...")
        
        # Get workflow name from WDL filename
        workflow_name = self.config.workflow_path.stem
        
        # Find and remove old run directories
        import glob
        pattern = f"[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_*_{workflow_name}"
        run_dirs = glob.glob(str(self.config.repo_root / pattern))
        
        for run_dir in run_dirs:
            try:
                import shutil
                shutil.rmtree(run_dir)
            except Exception as e:
                print(f"Warning: Failed to remove {run_dir}: {e}")
        
        print("✓ Cleanup done")
        return True
    
    def _backup_docker_config(self) -> bool:
        """Backup Docker config file and remove to prevent credential issues."""
        docker_config = Path.home() / ".docker" / "config.json"
        docker_config_bak = Path.home() / ".docker" / "config.json.bak"
        
        if docker_config.exists():
            try:
                # Only backup if not already backed up
                if not docker_config_bak.exists():
                    import shutil
                    shutil.copy(docker_config, docker_config_bak)
                    print("Docker config backed up")
                
                # Remove the original to prevent Docker from trying to use
                # credential helpers that may not be authenticated
                docker_config.unlink()
                print("Docker config removed to prevent credential issues")
                return True
            except Exception as e:
                print(f"Warning: Failed to backup/remove Docker config: {e}")
        return True
    
    def _restore_docker_config(self) -> bool:
        """Restore Docker config file from backup."""
        docker_config = Path.home() / ".docker" / "config.json"
        docker_config_bak = Path.home() / ".docker" / "config.json.bak"
        
        if docker_config_bak.exists():
            try:
                import shutil
                shutil.move(docker_config_bak, docker_config)
                print("Docker config restored")
                return True
            except Exception as e:
                print(f"Warning: Failed to restore Docker config: {e}")
        return True
    
    def _write_miniwdl_config(self) -> Path:
        """Write miniwdl configuration file.
        
        Returns:
            Path to written config file
        """
        self.config.tmp_dir.mkdir(parents=True, exist_ok=True)
        self.config.call_cache_dir.mkdir(parents=True, exist_ok=True)
        
        config_path = self.config.tmp_dir / "miniwdl.cfg"
        
        with open(config_path, "w") as f:
            f.write("[call_cache]\n")
            f.write("get = true\n")
            f.write("put = true\n")
            f.write(f"dir = {self.config.call_cache_dir}\n")
        
        return config_path
    
    def _setup_gcs_server(self) -> bool:
        """Check and setup fake-GCS server.
        
        Returns:
            True if GCS server is available and ready, False otherwise
        """
        print("Checking for fake-gcs-server...")
        
        if not self.gcs_manager.is_available():
            print("✗ fake-gcs-server is not running!")
            print("")
            print("Start it with:")
            print("  docker run -d -p 4443:4443 --name fake-gcs \\")
            print("    fsouza/fake-gcs-server -scheme http \\")
            print("    -external-url http://host.docker.internal:4443 \\")
            print("    -public-host host.docker.internal:4443")
            return False
        
        print("✓ fake-gcs-server is running")
        print("")
        
        # Ensure test data exists
        if not self.gcs_manager.ensure_test_data(self.config.reset_fake_gcs):
            return False
        
        print("")
        return True
    
    def _setup_and_shadow_docker(self) -> bool:
        """Build and shadow Docker image.
        
        Returns:
            True if successful, False otherwise
        """
        # Build test image (unless skipped)
        if not self.config.skip_build:
            if not self.docker_manager.build_test_image():
                return False
        else:
            print(f"Skipping Docker image build (SKIP_BUILD=true)")
            print(f"Using existing image {self.config.docker.test_image}")
        
        print("")
        
        # Shadow production image with test image
        if not self.docker_manager.shadow_prod_image():
            return False
        
        print("")
        return True
    
    def run(self) -> Tuple[bool, Optional[str]]:
        """Execute the complete test workflow.
        
        Returns:
            Tuple of (success: bool, error_message: Optional[str])
            - If successful: (True, None)
            - If failed: (False, error_message or path to stderr)
        """
        os.chdir(self.config.repo_root)
        
        print("=" * 70)
        print("WARP Test Harness - Workflow Execution")
        print("=" * 70)
        print("")
        
        # Log configuration
        print("Configuration:")
        for key, value in self.config.to_dict().items():
            if key not in ["docker", "gcs"]:
                print(f"  {key}: {value}")
        print("")
        
        try:
            # Cleanup old runs
            if not self.cleanup_old_runs():
                return False, "Cleanup failed"
            
            print("")
            
            # Setup GCS
            if not self._setup_gcs_server():
                return False, "fake-gcs-server unavailable"
            
            # Build and shadow Docker image BEFORE touching docker config.
            # On macOS, ~/.docker/config.json contains the currentContext needed
            # for Docker CLI to locate the Docker Desktop socket. Removing it
            # earlier causes all docker commands to fail.
            if not self._setup_and_shadow_docker():
                return False, "Docker setup failed"
            
            # Patch WDL for testing
            print("Patching WDL docker images for testing...")
            patched_result = self._patch_wdl()
            if not patched_result:
                self.docker_manager.restore_prod_image()
                return False, "WDL patching failed"
            
            patched_workflow = patched_result
            print("")
            
            # Write miniwdl config
            miniwdl_cfg = self._write_miniwdl_config()
            
            # Remove Docker config just before running miniwdl to prevent
            # credential helpers (e.g. docker-credential-gcloud) from being
            # invoked when Docker Swarm tries to authenticate against GCR.
            self._backup_docker_config()
            
            # Run workflow
            print(f"Running miniwdl workflow...")
            print(f"  Workflow: {self.config.workflow_wdl}")
            print(f"  Inputs: {self.config.inputs_json}")
            print("")
            
            success, error = self._run_miniwdl(patched_workflow, miniwdl_cfg)
            
            # Cleanup
            self.docker_manager.restore_prod_image()
            self._restore_docker_config()
            
            return success, error
        
        except Exception as e:
            print(f"✗ Unexpected error: {e}")
            self._restore_docker_config()
            self.docker_manager.restore_prod_image()
            return False, str(e)
    
    def _patch_wdl(self) -> Optional[Path]:
        """Patch WDL file for testing.
        
        Returns:
            Path to patched WDL if successful, None otherwise
        """
        replacements = create_simple_docker_replacements(
            self.config.docker.prod_image,
            self.config.docker.test_image,
        )
        
        patched_path = self.wdl_patcher.get_patched_path(self.config.workflow_path)
        
        if self.wdl_patcher.patch_workflow(
            self.config.workflow_path,
            patched_path,
            replacements,
        ):
            print(f"✓ WDL patched: {patched_path}")
            print("")
            return patched_path
        
        return None
    
    def _run_miniwdl(self, workflow_path: Path, config_path: Path) -> Tuple[bool, Optional[str]]:
        """Run miniwdl on the patched workflow.
        
        Args:
            workflow_path: Path to patched WDL file
            config_path: Path to miniwdl config file
        
        Returns:
            Tuple of (success: bool, error_message: Optional[str])
        """
        env = os.environ.copy()
        env["STORAGE_EMULATOR_HOST"] = self.config.gcs.storage_emulator_host
        
        # Set DOCKER_HOST based on configuration:
        # 1. If explicitly configured in config, use that
        # 2. If socket exists at ~/.docker/run/docker.sock (Linux), use that
        # 3. Otherwise let Docker CLI auto-detect (macOS Docker Desktop handles this)
        if self.config.docker.docker_host:
            env["DOCKER_HOST"] = self.config.docker.docker_host
        else:
            docker_socket = Path.home() / ".docker" / "run" / "docker.sock"
            if docker_socket.exists():
                env["DOCKER_HOST"] = f"unix://{docker_socket}"
        
        try:
            result = subprocess.run(
                [
                    "miniwdl",
                    "run",
                    str(workflow_path),
                    "-i", str(self.config.inputs_path),
                    "--cfg", str(config_path),
                    "--env", f"STORAGE_EMULATOR_HOST={self.config.gcs.storage_emulator_host}",
                ],
                env=env,
                check=False,
            )
            
            if result.returncode == 0:
                print("")
                print("✓ Workflow completed successfully!")
                return True, None
            else:
                print("")
                print(f"✗ Workflow failed with exit code {result.returncode}")
                
                # Try to find and show error from latest run
                error_output = self._collect_error_output()
                if error_output:
                    print("")
                    print("=" * 70)
                    print("ERROR OUTPUT")
                    print("=" * 70)
                    print(error_output)
                
                return False, error_output
        
        except FileNotFoundError:
            return False, "miniwdl command not found. Ensure miniwdl is installed."
    
    def _collect_error_output(self) -> Optional[str]:
        """Collect error output from latest failed miniwdl run.
        
        Returns:
            Error output if found, None otherwise
        """
        import glob
        
        workflow_name = self.config.workflow_path.stem
        pattern = f"{self.config.repo_root}/[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_*_{workflow_name}"
        
        # Find latest run directory
        runs = sorted(glob.glob(pattern), reverse=True)
        if not runs:
            return None
        
        latest_run = runs[0]
        
        # Find stderr files from failed tasks
        stderr_files = sorted(
            glob.glob(f"{latest_run}/**/stderr.txt", recursive=True),
            key=lambda f: os.path.getmtime(f),
            reverse=True,
        )
        
        if stderr_files:
            try:
                with open(stderr_files[0], "r") as f:
                    return f.read()
            except Exception:
                pass
        
        return None
