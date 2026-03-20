"""
Command-line interface for WARP test harness.

Provides CLI for running tests with various options and configurations.
"""

import argparse
import sys
from pathlib import Path
from . import __version__
from .config import get_pipeline_config, create_custom_config
from .runner import WorkflowRunner


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser for CLI.
    
    Returns:
        ArgumentParser configured for test harness
    """
    parser = argparse.ArgumentParser(
        prog="warp-test-harness",
        description="WARP workflow test harness - Run WDL pipelines with Docker and fake-GCS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run predefined mt_coverage_merge tests
  python -m scripts.test_harness mt_coverage_merge

  # Run tests with options
  python -m scripts.test_harness mt_coverage_merge --skip-cleanup --reset-gcs
  
  # Run custom pipeline
  python -m scripts.test_harness custom \\
    --workflow pipelines/wdl/my_pipeline/pipeline.wdl \\
    --inputs testing/my_inputs.json \\
    --prod-image "us.gcr.io/my-org/my-image:v1.0" \\
    --test-image "my-image:local-test" \\
    --dockerfile testing/my-image/

Notes:
  - fake-gcs-server must be running before test execution
  - miniwdl and Docker must be installed
  - Run from repository root or use --repo-root to specify location
        """,
    )
    
    parser.add_argument(
        "pipeline",
        nargs="?",
        help=(
            "Pipeline to test (e.g., 'mt_coverage_merge') or 'custom' for custom configuration. "
            "Use --list to see available pipelines."
        ),
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available predefined pipelines",
    )
    
    parser.add_argument(
        "--repo-root",
        type=str,
        default=".",
        help="Repository root directory (default: current directory)",
    )
    
    # Custom pipeline options
    parser.add_argument(
        "--workflow",
        type=str,
        help="Path to WDL workflow (relative to repo root) - required for 'custom' pipeline",
    )
    
    parser.add_argument(
        "--inputs",
        type=str,
        help="Path to inputs JSON (relative to repo root) - required for 'custom' pipeline",
    )
    
    parser.add_argument(
        "--prod-image",
        type=str,
        help="Production Docker image - required for 'custom' pipeline",
    )
    
    parser.add_argument(
        "--test-image",
        type=str,
        help="Test Docker image - required for 'custom' pipeline",
    )
    
    parser.add_argument(
        "--dockerfile",
        type=str,
        help="Path to Dockerfile directory (relative to repo root) - optional for 'custom' pipeline",
    )
    
    parser.add_argument(
        "--docker-host",
        type=str,
        help="Docker host connection string (e.g., 'unix:///var/run/docker.sock'). Auto-detected if not provided.",
    )
    
    # Control flags
    parser.add_argument(
        "--skip-cleanup",
        action="store_true",
        help="Skip cleanup of old run directories (useful for inspecting outputs)",
    )
    
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip Docker image build (use existing image)",
    )
    
    parser.add_argument(
        "--reset-gcs",
        action="store_true",
        help="Force re-setup of fake-GCS data (overwrites existing)",
    )
    
    return parser


def list_pipelines() -> None:
    """List available predefined pipelines."""
    from .config import PIPELINE_CONFIGS
    
    print("Available pipelines:")
    for name in sorted(PIPELINE_CONFIGS.keys()):
        print(f"  - {name}")
    print("")
    print("Use 'custom' for custom pipeline configuration.")


def main(argv: list = None) -> int:
    """Main entry point for CLI.
    
    Args:
        argv: Command-line arguments (defaults to sys.argv[1:])
    
    Returns:
        Exit code (0 for success, 1 for failure)
    """
    parser = create_parser()
    args = parser.parse_args(argv)
    
    # Handle --list flag
    if args.list:
        list_pipelines()
        return 0
    
    # Check if pipeline is provided
    if not args.pipeline:
        parser.print_help()
        return 1
    
    # Get configuration
    try:
        if args.pipeline == "custom":
            # Validate required custom pipeline arguments
            required = ["workflow", "inputs", "prod_image", "test_image"]
            missing = [arg for arg in required if not getattr(args, arg.lower())]
            if missing:
                missing_flags = ", ".join(f"--{arg.replace('_', '-')}" for arg in missing)
                parser.error(f"Custom pipeline requires: {missing_flags}")
            
            config = create_custom_config(
                repo_root=args.repo_root,
                workflow_wdl=args.workflow,
                inputs_json=args.inputs,
                prod_image=args.prod_image,
                test_image=args.test_image,
                dockerfile_dir=args.dockerfile,
                docker_host=args.docker_host,
                skip_cleanup=args.skip_cleanup,
                skip_build=args.skip_build,
                reset_fake_gcs=args.reset_gcs,
            )
        else:
            config = get_pipeline_config(args.pipeline, args.repo_root)
            config.skip_cleanup = args.skip_cleanup
            config.skip_build = args.skip_build
            config.reset_fake_gcs = args.reset_gcs
            if args.docker_host:
                config.docker.docker_host = args.docker_host
    
    except ValueError as e:
        parser.error(str(e))
        return 1
    except FileNotFoundError as e:
        parser.error(str(e))
        return 1
    
    # Run workflow
    runner = WorkflowRunner(config)
    success, error = runner.run()
    
    if success:
        return 0
    else:
        if error:
            print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
