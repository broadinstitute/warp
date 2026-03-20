"""
Configuration management for WARP test harness.

Defines test configuration with support for different pipelines and environments.
All paths are resolved relative to the repository root.
"""

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, List


@dataclass
class DockerConfig:
    """Docker image configuration for test execution."""
    
    # Production image that will be shadowed during testing
    prod_image: str
    
    # Local test image (usually built for testing)
    test_image: str
    
    # Optional: dockerfile path relative to repo root (for building test image)
    dockerfile_dir: Optional[str] = None
    
    # Optional: Docker host connection string (e.g., "unix:///var/run/docker.sock")
    # If None, Docker CLI will auto-detect (recommended for macOS Docker Desktop)
    docker_host: Optional[str] = None


@dataclass
class GCSConfig:
    """Fake-GCS server configuration."""
    
    # Fake-GCS server URL
    host: str = "http://localhost:4443"
    
    # Bucket name for test data
    bucket: str = "fake-bucket"
    
    # Storage emulator host (for environment variable in containers)
    storage_emulator_host: str = "http://host.docker.internal:4443"
    
    # Test objects required to exist in fake-GCS
    required_objects: List[str] = field(default_factory=lambda: [
        "coverage_s001.tsv",
        "coverage_s002.tsv",
        "s001.vcf",
        "s002.vcf",
        "vcf_shard_00000_shard_mt.tar.gz",
        "blacklist_sites.hg38.chrM.bed",
    ])


@dataclass
class TestConfig:
    """Test configuration for a WARP pipeline."""
    
    # Paths (relative to REPO_ROOT)
    workflow_wdl: str
    inputs_json: str
    
    # Docker configuration
    docker: DockerConfig
    
    # GCS configuration
    gcs: GCSConfig = field(default_factory=GCSConfig)
    
    # Repository root (absolute path)
    repo_root: Path = field(default_factory=Path.cwd)
    
    # Temporary directory for test artifacts
    tmp_dir: Optional[Path] = None
    
    # Call cache directory
    call_cache_dir: Optional[Path] = None
    
    # Control flags
    skip_cleanup: bool = False
    skip_build: bool = False
    reset_fake_gcs: bool = False
    
    def __post_init__(self):
        """Resolve all paths to absolute paths."""
        # Convert repo_root to Path if string
        if isinstance(self.repo_root, str):
            self.repo_root = Path(self.repo_root)
        
        # Set defaults for tmp and cache directories
        if self.tmp_dir is None:
            self.tmp_dir = self.repo_root / "tmp"
        elif isinstance(self.tmp_dir, str):
            self.tmp_dir = Path(self.tmp_dir)
        
        if self.call_cache_dir is None:
            self.call_cache_dir = self.tmp_dir / "miniwdl-call-cache"
        elif isinstance(self.call_cache_dir, str):
            self.call_cache_dir = Path(self.call_cache_dir)
    
    @property
    def workflow_path(self) -> Path:
        """Absolute path to workflow WDL."""
        wdl_path = self.repo_root / self.workflow_wdl
        if not wdl_path.exists():
            raise FileNotFoundError(f"Workflow not found: {wdl_path}")
        return wdl_path
    
    @property
    def inputs_path(self) -> Path:
        """Absolute path to inputs JSON."""
        inputs_path = self.repo_root / self.inputs_json
        if not inputs_path.exists():
            raise FileNotFoundError(f"Inputs JSON not found: {inputs_path}")
        return inputs_path
    
    @property
    def dockerfile_path(self) -> Optional[Path]:
        """Absolute path to Dockerfile directory (for building test image)."""
        if self.docker.dockerfile_dir is None:
            return None
        dockerfile_path = self.repo_root / self.docker.dockerfile_dir
        if not dockerfile_path.exists():
            raise FileNotFoundError(f"Dockerfile directory not found: {dockerfile_path}")
        return dockerfile_path
    
    def to_dict(self) -> Dict:
        """Convert config to dictionary for logging/debugging."""
        return {
            "workflow_wdl": self.workflow_wdl,
            "inputs_json": self.inputs_json,
            "docker": {
                "prod_image": self.docker.prod_image,
                "test_image": self.docker.test_image,
                "dockerfile_dir": self.docker.dockerfile_dir,
                "docker_host": self.docker.docker_host,
            },
            "gcs": {
                "host": self.gcs.host,
                "bucket": self.gcs.bucket,
                "storage_emulator_host": self.gcs.storage_emulator_host,
                "required_objects": self.gcs.required_objects,
            },
            "repo_root": str(self.repo_root),
            "tmp_dir": str(self.tmp_dir),
            "call_cache_dir": str(self.call_cache_dir),
            "skip_cleanup": self.skip_cleanup,
            "skip_build": self.skip_build,
            "reset_fake_gcs": self.reset_fake_gcs,
        }


# Predefined configurations for known pipelines
PIPELINE_CONFIGS = {
    "mt_coverage_merge": lambda repo_root: TestConfig(
        workflow_wdl="all_of_us/mitochondria/mt_coverage_merge.wdl",
        inputs_json="all_of_us/mitochondria/testing/inputs/mt_coverage_merge_inputs.json",
        docker=DockerConfig(
            prod_image="us.gcr.io/broad-gotc-prod/aou-mitochondrial-combine-vcfs-covdb:1.0.1",
            test_image="aou-mitochondrial-combine-vcfs-covdb:local-test",
            dockerfile_dir="all_of_us/mitochondria/testing/dockers/aou-mitochondrial-combine-vcfs-covdb",
        ),
        gcs=GCSConfig(),
        repo_root=repo_root,
    ),
}


def get_pipeline_config(pipeline_name: str, repo_root: str = ".") -> TestConfig:
    """Get test configuration for a named pipeline.
    
    Args:
        pipeline_name: Name of the pipeline (e.g., 'mt_coverage_merge')
        repo_root: Repository root directory (default: current directory)
    
    Returns:
        TestConfig configured for the specified pipeline
    
    Raises:
        ValueError: If pipeline_name is not found in PIPELINE_CONFIGS
    """
    if pipeline_name not in PIPELINE_CONFIGS:
        available = ", ".join(PIPELINE_CONFIGS.keys())
        raise ValueError(
            f"Unknown pipeline: {pipeline_name}. Available: {available}"
        )
    
    if repo_root == ".":
        repo_root = os.getcwd()
    
    return PIPELINE_CONFIGS[pipeline_name](repo_root)


def create_custom_config(
    repo_root: str,
    workflow_wdl: str,
    inputs_json: str,
    prod_image: str,
    test_image: str,
    dockerfile_dir: Optional[str] = None,
    docker_host: Optional[str] = None,
    **kwargs
) -> TestConfig:
    """Create a custom test configuration for any pipeline.
    
    Args:
        repo_root: Repository root directory
        workflow_wdl: Path to WDL workflow (relative to repo root)
        inputs_json: Path to inputs JSON (relative to repo root)
        prod_image: Production Docker image
        test_image: Test Docker image
        dockerfile_dir: Optional path to Dockerfile directory (relative to repo root)
        docker_host: Optional Docker host connection string
        **kwargs: Additional TestConfig parameters (skip_cleanup, etc.)
    
    Returns:
        Custom TestConfig instance
    """
    return TestConfig(
        workflow_wdl=workflow_wdl,
        inputs_json=inputs_json,
        docker=DockerConfig(
            prod_image=prod_image,
            test_image=test_image,
            dockerfile_dir=dockerfile_dir,
            docker_host=docker_host,
        ),
        repo_root=repo_root,
        **kwargs
    )
