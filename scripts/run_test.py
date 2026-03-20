#!/usr/bin/env python3
"""
WARP Test Harness - Entry point script.

Run WARP WDL pipelines with Docker and fake-GCS server support.

This script is designed to run on the HOST MACHINE (macOS/Linux with Docker and miniwdl)
and NOT from within the dev container.

Examples:
  # Run mt_coverage_merge tests
  ./scripts/run_test.py mt_coverage_merge

  # Run with custom options
  ./scripts/run_test.py mt_coverage_merge --skip-cleanup --reset-gcs

  # Run custom pipeline
  ./scripts/run_test.py custom \\
    --workflow pipelines/wdl/my_pipeline/pipeline.wdl \\
    --inputs testing/my_inputs.json \\
    --prod-image "us.gcr.io/my-org/my-image:v1.0" \\
    --test-image "my-image:local-test"

  # List available pipelines
  ./scripts/run_test.py --list

Environment Variables:
  WARP_SKIP_CLEANUP    Set to "true" to preserve old run directories
  WARP_SKIP_BUILD      Set to "true" to skip Docker image build
  WARP_RESET_GCS       Set to "true" to force reset of fake-GCS data
"""

import sys
import os
from pathlib import Path

# Add parent directory to path so we can import test_harness
scripts_dir = Path(__file__).parent
sys.path.insert(0, str(scripts_dir))

from test_harness.cli import main


def main_with_env() -> int:
    """Main entry point with environment variable support.
    
    Converts environment variables to command-line arguments.
    """
    argv = sys.argv[1:]
    
    # Convert environment variables to command-line flags
    if os.getenv("WARP_SKIP_CLEANUP", "").lower() == "true":
        argv.append("--skip-cleanup")
    
    if os.getenv("WARP_SKIP_BUILD", "").lower() == "true":
        argv.append("--skip-build")
    
    if os.getenv("WARP_RESET_GCS", "").lower() == "true":
        argv.append("--reset-gcs")
    
    return main(argv)


if __name__ == "__main__":
    sys.exit(main_with_env())
