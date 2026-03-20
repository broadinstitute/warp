"""
Docker image management for test harness.

Handles building, shadowing, and restoring Docker images for testing.
"""

import subprocess
import json
import os
from pathlib import Path
from typing import Optional
from .config import DockerConfig


class DockerManager:
    """Manages Docker image lifecycle for testing."""
    
    def __init__(self, config: DockerConfig):
        """Initialize Docker manager.
        
        Args:
            config: DockerConfig with image names and paths
        """
        self.config = config
        self.original_image_id: Optional[str] = None
    
    def _get_docker_env(self) -> dict:
        """Get environment dict with DOCKER_HOST configured if needed.
        
        NOTE: This method is kept for reference/future use but is NOT currently used
        for docker build/tag/rmi commands. Those commands are allowed to inherit the
        parent environment to benefit from Docker Desktop's auto-detection on macOS
        and similar behaviors on other platforms.
        
        Returns:
            Modified environment dict for docker commands
        """
        env = os.environ.copy()
        
        # Set DOCKER_HOST if configured
        if self.config.docker_host:
            env["DOCKER_HOST"] = self.config.docker_host
        else:
            # Try to detect Docker Desktop socket on macOS/Linux
            docker_socket = Path.home() / ".docker" / "run" / "docker.sock"
            if docker_socket.exists():
                env["DOCKER_HOST"] = f"unix://{docker_socket}"
        
        return env
    
    def build_test_image(self) -> bool:
        """Build local test Docker image.
        
        Returns:
            True if build succeeded, False otherwise
        """
        if self.config.dockerfile_dir is None:
            print(f"Dockerfile directory not configured, skipping build")
            return True
        
        print(f"Building local test image ({self.config.test_image})...")
        try:
            # Don't override environment - let docker CLI use whatever works in parent shell
            # (Docker Desktop auto-detection, etc.)
            result = subprocess.run(
                ["docker", "build", "-t", self.config.test_image, str(self.config.dockerfile_dir)],
                capture_output=True,
                text=True,
                check=False,
            )
            
            if result.returncode != 0:
                print(f"Docker build failed:")
                print(result.stderr)
                return False
            
            # Print last few lines of output
            lines = result.stdout.strip().split("\n")
            for line in lines[-5:]:
                if line.strip():
                    print(line)
            
            print(f"✓ Test image built successfully")
            return True
        except FileNotFoundError:
            print("✗ Docker command not found. Ensure Docker is installed and on PATH.")
            return False
    
    def get_image_id(self, image_name: str) -> Optional[str]:
        """Get Docker image ID for the given image name.
        
        Args:
            image_name: Docker image name or tag
        
        Returns:
            Image ID if found, None otherwise
        """
        try:
            result = subprocess.run(
                ["docker", "inspect", "--format={{.Id}}", image_name],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode == 0:
                return result.stdout.strip()
        except FileNotFoundError:
            pass
        return None
    
    def shadow_prod_image(self) -> bool:
        """Replace production image with test image.
        
        Removes the production image tag and retags the test image with it.
        This ensures miniwdl uses the test image regardless of lookup method.
        
        Returns:
            True if successful, False otherwise
        """
        print(f"Shadowing production image with test image...")
        
        # Save the original image ID if it exists
        self.original_image_id = self.get_image_id(self.config.prod_image)
        
        if self.original_image_id:
            print(f"Removing old {self.config.prod_image}...")
            try:
                subprocess.run(
                    ["docker", "rmi", self.config.prod_image],
                    capture_output=True,
                    check=False,
                )
            except FileNotFoundError:
                return False
        
        # Tag test image as production image
        print(f"Tagging test image as {self.config.prod_image}...")
        try:
            result = subprocess.run(
                ["docker", "tag", self.config.test_image, self.config.prod_image],
                capture_output=True,
                text=True,
                check=False,
            )
            
            if result.returncode != 0:
                print(f"Failed to tag image: {result.stderr}")
                return False
            
            print(f"✓ Production image replaced with test image")
            return True
        except FileNotFoundError:
            return False
    
    def restore_prod_image(self) -> bool:
        """Restore production image tag.
        
        Removes the test image tag (production name) so it doesn't shadow.
        Note: The original image is not restored (it was removed); user may need
        to re-pull it on next run if needed.
        
        Returns:
            True if successful, False otherwise
        """
        print(f"Cleaning up image shadowing...")
        
        try:
            # Remove the test image tag
            subprocess.run(
                ["docker", "rmi", self.config.prod_image],
                capture_output=True,
                check=False,
            )
            
            if self.original_image_id:
                print(f"Original image was replaced; you may need to re-pull on next run")
            
            print(f"✓ Production image tag cleaned up")
            return True
        except FileNotFoundError:
            return False
