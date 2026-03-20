"""
Fake-GCS server management for test harness.

Handles validation and setup of test data in fake-GCS server.
"""

import subprocess
import time
from pathlib import Path
from typing import List, Optional
from .config import GCSConfig


class GCSManager:
    """Manages fake-GCS server setup and test data."""
    
    def __init__(self, config: GCSConfig, repo_root: Path):
        """Initialize GCS manager.
        
        Args:
            config: GCSConfig with server settings
            repo_root: Repository root path (for locating setup scripts)
        """
        self.config = config
        self.repo_root = Path(repo_root)
        self.setup_script = self.repo_root / "testing" / "setup_fake_gcs.sh"
    
    def is_available(self, timeout: float = 5.0) -> bool:
        """Check if fake-GCS server is running and available.
        
        Args:
            timeout: Timeout in seconds for the check
        
        Returns:
            True if server is available, False otherwise
        """
        try:
            result = subprocess.run(
                ["curl", "-s", "-f", "-o", "/dev/null", f"{self.config.host}/storage/v1/b"],
                timeout=timeout,
                capture_output=True,
                check=False,
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False
    
    def check_required_objects(self) -> bool:
        """Check if required test objects exist in fake-GCS.
        
        Returns:
            True if all required objects exist, False if any are missing
        """
        for obj in self.config.required_objects:
            if not self._object_exists(obj):
                return False
        return True
    
    def _object_exists(self, object_name: str) -> bool:
        """Check if a specific object exists in fake-GCS bucket.
        
        Args:
            object_name: Name of the object to check
        
        Returns:
            True if object exists, False otherwise
        """
        try:
            url = f"{self.config.host}/storage/v1/b/{self.config.bucket}/o/{object_name}"
            result = subprocess.run(
                ["curl", "-s", "-f", "-o", "/dev/null", url],
                timeout=2.0,
                capture_output=True,
                check=False,
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False
    
    def setup_test_data(self) -> bool:
        """Run fake-GCS setup script to upload test data.
        
        Returns:
            True if setup succeeded, False otherwise
        """
        if not self.setup_script.exists():
            print(f"Setup script not found: {self.setup_script}")
            return False
        
        print(f"Running fake-GCS setup to upload test data...")
        try:
            result = subprocess.run(
                ["bash", str(self.setup_script), str(self.repo_root)],
                capture_output=True,
                text=True,
                check=False,
            )
            
            if result.returncode != 0:
                print(f"Setup script failed:")
                print(result.stderr)
                return False
            
            # Print output
            if result.stdout:
                for line in result.stdout.strip().split("\n"):
                    if line.strip():
                        print(line)
            
            print(f"✓ Test data setup complete")
            return True
        except FileNotFoundError:
            print("bash command not found")
            return False
    
    def ensure_test_data(self, force_reset: bool = False) -> bool:
        """Ensure test data exists in fake-GCS, uploading if necessary.
        
        Args:
            force_reset: If True, always re-run setup (overwrites existing data)
        
        Returns:
            True if data is available, False if setup failed
        """
        if force_reset:
            print("Force resetting fake-GCS data (RESET_FAKE_GCS=true)...")
            return self.setup_test_data()
        
        # Check if required objects already exist
        missing_objects = [
            obj for obj in self.config.required_objects
            if not self._object_exists(obj)
        ]
        
        if not missing_objects:
            print(f"✓ Test data already available in fake-GCS")
            return True
        
        print(f"Test data missing from fake-GCS (missing: {', '.join(missing_objects[:3])}...)")
        return self.setup_test_data()
