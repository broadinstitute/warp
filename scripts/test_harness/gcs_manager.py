"""
Fake-GCS server management for test harness.

Handles validation and setup of test data in fake-GCS server.
"""

import subprocess
import tarfile
import tempfile
import urllib.parse
from pathlib import Path
from typing import List, Optional, Tuple
from .config import GCSConfig


class GCSManager:
    """Manages fake-GCS server setup and test data."""
    
    def __init__(self, config: GCSConfig, repo_root: Path):
        """Initialize GCS manager.
        
        Args:
            config: GCSConfig with server settings
            repo_root: Repository root path (for locating test data files)
        """
        self.config = config
        self.repo_root = Path(repo_root)
    
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
    
    def _create_bucket(self) -> bool:
        """Create the test bucket in fake-GCS.
        
        Returns:
            True if bucket created or already exists, False on error
        """
        try:
            result = subprocess.run(
                [
                    "curl",
                    "-s",
                    "-X",
                    "POST",
                    f"{self.config.host}/storage/v1/b",
                    "-H",
                    "Content-Type: application/json",
                    "-d",
                    f'{{"name":"{self.config.bucket}"}}',
                ],
                capture_output=True,
                text=True,
                check=False,
                timeout=5.0,
            )
            if result.returncode == 0:
                print("✓ Bucket created (or already exists)")
                return True
            else:
                print("⚠ Could not create bucket (may already exist)")
                return True
        except subprocess.TimeoutExpired:
            print("✗ Bucket creation timed out")
            return False
        except FileNotFoundError:
            print("✗ curl command not found")
            return False

    def _upload_file(self, local_path: Path, gcs_object_name: str) -> bool:
        """Upload a file to fake-GCS bucket.
        
        Args:
            local_path: Path to local file
            gcs_object_name: GCS object name (can include paths like "resources/dir/file.txt")
        
        Returns:
            True if upload succeeded, False otherwise
        """
        if not local_path.exists():
            print(f"✗ Missing file: {local_path}")
            return False

        try:
            # URL encode the object name (preserve slashes for directories)
            encoded_name = urllib.parse.quote(gcs_object_name, safe="")
            
            url = f"{self.config.host}/upload/storage/v1/b/{self.config.bucket}/o?uploadType=media&name={encoded_name}"
            
            with open(local_path, "rb") as f:
                result = subprocess.run(
                    ["curl", "-s", "-X", "POST", url, "--data-binary", "@-"],
                    stdin=f,
                    capture_output=True,
                    check=False,
                    timeout=10.0,
                )
            
            if result.returncode == 0:
                print(f"  ✓ Uploaded {gcs_object_name}")
                return True
            else:
                print(f"  ✗ Failed to upload {gcs_object_name}")
                return False
        except subprocess.TimeoutExpired:
            print(f"  ✗ Upload timed out: {gcs_object_name}")
            return False
        except (FileNotFoundError, OSError) as e:
            print(f"  ✗ Upload error: {e}")
            return False

    def _create_mt_shard_tar(self) -> Tuple[bool, Optional[Path]]:
        """Create MT shard tar file for merge tests.
        
        Returns:
            Tuple of (success: bool, tar_path: Optional[Path])
        """
        try:
            # Create temporary structure
            tmp_dir = Path(tempfile.gettempdir()) / "warp_mt_shards"
            tmp_dir.mkdir(parents=True, exist_ok=True)
            
            shard_dir = tmp_dir / "vcf_shard_00000"
            shard_dir.mkdir(exist_ok=True)
            
            # Create minimal MT VCF file
            vcf_file = shard_dir / "vcf_shard_00000.mt"
            vcf_file.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "MT\t100\t.\tA\tG\t.\tPASS\t.\n"
            )
            
            # Create tar.gz
            tar_path = tmp_dir / "vcf_shard_00000_shard_mt.tar.gz"
            with tarfile.open(tar_path, "w:gz") as tar:
                tar.add(shard_dir, arcname="vcf_shard_00000")
            
            print("  ✓ Created vcf_shard_00000_shard_mt.tar.gz")
            return True, tar_path
        except Exception as e:
            print(f"  ✗ Failed to create MT shard tar: {e}")
            return False, None

    def setup_test_data(self) -> bool:
        """Setup fake-GCS test data by uploading all required files.
        
        This reproduces all features from all_of_us/mitochondria/testing/setup_fake_gcs.sh in Python:
        1. Create bucket
        2. Upload coverage TSV files
        3. Upload VCF files
        4. Create and upload MT shard tar files
        5. Upload BED files
        6. Upload mock resource files with proper GCS paths
        
        Returns:
            True if setup succeeded, False otherwise
        """
        print("Setting up fake-GCS server for testing...")
        print(f"Repository root: {self.repo_root}")
        print()
        
        # Step 1: Create bucket
        print("Creating bucket...", end=" ")
        if not self._create_bucket():
            return False
        print()

        # Step 2: Upload coverage files (coverage_s*.tsv)
        print("Uploading coverage files...")
        coverage_files = [
            self.repo_root / "all_of_us" / "mitochondria" / "testing" / "mocks" / "coverage_s001.tsv",
            self.repo_root / "all_of_us" / "mitochondria" / "testing" / "mocks" / "coverage_s002.tsv",
        ]
        for coverage_file in coverage_files:
            if not self._upload_file(coverage_file, coverage_file.name):
                return False
        print()

        # Step 3: Upload VCF files (s*.vcf)
        print("Uploading VCF files...")
        vcf_files = [
            self.repo_root / "all_of_us" / "mitochondria" / "testing" / "mocks" / "s001.vcf",
            self.repo_root / "all_of_us" / "mitochondria" / "testing" / "mocks" / "s002.vcf",
        ]
        for vcf_file in vcf_files:
            if not self._upload_file(vcf_file, vcf_file.name):
                return False
        print()

        # Step 4: Create and upload MT shard tar files
        print("Creating MT shard tar files...")
        success, tar_path = self._create_mt_shard_tar()
        if not success or not tar_path:
            return False
        if not self._upload_file(tar_path, tar_path.name):
            return False
        print()

        # Step 5: Upload BED file (blacklist sites)
        print("Uploading artifact-prone-sites BED file...")
        blacklist_file = self.repo_root / "all_of_us" / "mitochondria" / "testing" / "mocks" / "blacklist_sites.hg38.chrM.bed"
        if not self._upload_file(blacklist_file, "blacklist_sites.hg38.chrM.bed"):
            return False
        print()

        # Step 6: Upload mock resource files with specific GCS paths
        print("Uploading mock resource files for annotations...")
        mock_resources = [
            (
                self.repo_root / "all_of_us" / "mitochondria" / "testing" / "mocks" / "phylotree_mock.tsv",
                "resources/mitochondria/phylotree/rCRS-centered_phylo_vars_final_update.txt",
            ),
            (
                self.repo_root / "all_of_us" / "mitochondria" / "testing" / "mocks" / "variant_context_mock.tsv",
                "resources/mitochondria/variant_context/chrM_pos_ref_alt_context_categories.txt",
            ),
            (
                self.repo_root / "all_of_us" / "mitochondria" / "testing" / "mocks" / "pon_mt_trna_mock.tsv",
                "resources/mitochondria/trna_predictions/pon_mt_trna_predictions_08_27_2020.txt",
            ),
            (
                self.repo_root / "all_of_us" / "mitochondria" / "testing" / "mocks" / "mitotip_mock.tsv",
                "resources/mitochondria/trna_predictions/mitotip_scores_08_27_2020.txt",
            ),
        ]
        for local_path, gcs_path in mock_resources:
            if not self._upload_file(local_path, gcs_path):
                return False
        print()

        print("✓ Setup complete!")
        return True
    
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
