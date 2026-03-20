"""
WDL file patching for test harness.

Handles modifications to WDL files for testing (e.g., replacing docker images).
"""

import re
import shutil
from pathlib import Path
from typing import Optional


class WDLPatcher:
    """Patches WDL files for test execution."""
    
    def __init__(self, repo_root: Path):
        """Initialize WDL patcher.
        
        Args:
            repo_root: Repository root path
        """
        self.repo_root = Path(repo_root)
    
    def patch_workflow(
        self,
        workflow_path: Path,
        output_path: Path,
        image_replacements: Optional[dict] = None,
    ) -> bool:
        """Patch a WDL workflow file for testing.
        
        Args:
            workflow_path: Path to original WDL file
            output_path: Path to write patched WDL
            image_replacements: Dict mapping production image patterns to test images
                                e.g., {"prod-image:v1": "test-image:local-test"}
        
        Returns:
            True if patching succeeded, False otherwise
        """
        if not workflow_path.exists():
            print(f"✗ WDL file not found: {workflow_path}")
            return False
        
        try:
            # Ensure output directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Read original WDL
            with open(workflow_path, "r") as f:
                content = f.read()
            
            # Apply replacements
            if image_replacements:
                for prod_image, test_image in image_replacements.items():
                    # Replace docker image references
                    # E.g., "us.gcr.io/broad-gotc-prod/image:1.0.1" -> "image:local-test"
                    content = self._replace_docker_image(content, prod_image, test_image)
            
            # Write patched WDL
            with open(output_path, "w") as f:
                f.write(content)
            
            print(f"✓ WDL patched: {output_path}")
            return True
        except Exception as e:
            print(f"✗ Failed to patch WDL: {e}")
            return False
    
    def _replace_docker_image(
        self,
        content: str,
        prod_image: str,
        test_image: str,
    ) -> str:
        """Replace docker image reference in WDL content.
        
        Handles both fully-qualified image names and version tag replacements.
        
        Args:
            content: WDL file content
            prod_image: Production image pattern (may include version wildcards)
            test_image: Test image to replace with
        
        Returns:
            Modified WDL content
        """
        # If prod_image contains version like "image:v1.0", replace exact match
        if ":" in prod_image:
            # Exact version replacement (e.g., "us.gcr.io/.../:1.0.1" -> "test:local")
            content = content.replace(prod_image, test_image)
        else:
            # Pattern replacement - replace image with any version tag
            # Extract the image name and replace with version wildcard
            escaped_image = re.escape(prod_image)
            # Match image:any-version pattern
            pattern = f"{escaped_image}:[^\\s\"']*"
            content = re.sub(pattern, test_image, content)
        
        return content
    
    def get_patched_path(self, workflow_path: Path) -> Path:
        """Get the expected output path for a patched WDL file.
        
        Args:
            workflow_path: Path to original WDL file
        
        Returns:
            Path where patched file would be written (in tmp/)
        """
        tmp_dir = self.repo_root / "tmp"
        stem = workflow_path.stem
        return tmp_dir / f"{stem}.test.wdl"


def create_simple_docker_replacements(
    prod_image: str,
    test_image: str,
) -> dict:
    """Create a simple docker replacement mapping.
    
    Args:
        prod_image: Production image (full name or with version)
        test_image: Test image replacement
    
    Returns:
        Dict suitable for patch_workflow's image_replacements parameter
    """
    return {prod_image: test_image}
