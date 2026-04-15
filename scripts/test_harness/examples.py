"""
Examples of using the WARP test harness programmatically.

This file demonstrates various ways to configure and run tests
using the test harness Python API.
"""

from pathlib import Path
from scripts.test_harness.config import (
    TestConfig,
    DockerConfig,
    GCSConfig,
    get_pipeline_config,
    create_custom_config,
)
from scripts.test_harness.runner import WorkflowRunner


def example_1_predefined_pipeline():
    """Example 1: Run predefined pipeline configuration."""
    print("=" * 70)
    print("Example 1: Run predefined pipeline")
    print("=" * 70)
    
    # Get predefined configuration for mt_coverage_merge
    config = get_pipeline_config("mt_coverage_merge", repo_root=".")
    
    # Run the workflow
    runner = WorkflowRunner(config)
    success, error = runner.run()
    
    if success:
        print("✓ Test passed!")
    else:
        print(f"✗ Test failed: {error}")


def example_2_custom_pipeline():
    """Example 2: Run custom pipeline with manual configuration."""
    print("=" * 70)
    print("Example 2: Run custom pipeline")
    print("=" * 70)
    
    # Create custom configuration
    config = create_custom_config(
        repo_root=".",
        workflow_wdl="pipelines/wdl/my_pipeline/pipeline.wdl",
        inputs_json="testing/my_pipeline_inputs.json",
        prod_image="us.gcr.io/my-org/my-image:1.0.0",
        test_image="my-image:local-test",
        dockerfile_dir="testing/my-image/",
        skip_cleanup=True,  # Keep run directories for inspection
    )
    
    # Run the workflow
    runner = WorkflowRunner(config)
    success, error = runner.run()
    
    if success:
        print("✓ Test passed!")
    else:
        print(f"✗ Test failed: {error}")


def example_3_full_control():
    """Example 3: Full control with TestConfig and component managers."""
    print("=" * 70)
    print("Example 3: Full control with component managers")
    print("=" * 70)
    
    from scripts.test_harness.docker_manager import DockerManager
    from scripts.test_harness.gcs_manager import GCSManager
    from scripts.test_harness.wdl_patcher import WDLPatcher
    
    # Create detailed configuration
    config = TestConfig(
        workflow_wdl="pipelines/wdl/my_pipeline/pipeline.wdl",
        inputs_json="testing/my_pipeline_inputs.json",
        docker=DockerConfig(
            prod_image="us.gcr.io/my-org/my-image:1.0.0",
            test_image="my-image:local-test",
            dockerfile_dir="testing/my-image/",
        ),
        gcs=GCSConfig(
            host="http://localhost:4443",
            bucket="fake-bucket",
            storage_emulator_host="http://host.docker.internal:4443",
        ),
        repo_root=".",
        skip_cleanup=False,
        skip_build=False,
        reset_fake_gcs=False,
    )
    
    # Can also use individual managers
    docker_manager = DockerManager(config.docker)
    gcs_manager = GCSManager(config.gcs, config.repo_root)
    wdl_patcher = WDLPatcher(config.repo_root)
    
    print(f"Configuration:")
    print(f"  Workflow: {config.workflow_wdl}")
    print(f"  Inputs: {config.inputs_json}")
    print(f"  Prod image: {config.docker.prod_image}")
    print(f"  Test image: {config.docker.test_image}")
    print()
    
    # Check prerequisites
    print("Checking prerequisites...")
    if not gcs_manager.is_available():
        print("✗ fake-gcs-server is not running!")
        print("  Start it with:")
        print("    docker run -d -p 4443:4443 --name fake-gcs \\")
        print("      fsouza/fake-gcs-server -scheme http \\")
        print("      -external-url http://host.docker.internal:4443 \\")
        print("      -public-host host.docker.internal:4443")
        return
    
    print("✓ fake-gcs-server is running")
    print()
    
    # Ensure test data
    print("Ensuring test data in fake-GCS...")
    if gcs_manager.ensure_test_data():
        print("✓ Test data ready")
    else:
        print("✗ Failed to setup test data")
        return
    
    print()
    
    # Build Docker image
    print("Building test Docker image...")
    if docker_manager.build_test_image():
        print("✓ Docker image built")
    else:
        print("✗ Failed to build Docker image")
        return
    
    print()
    
    # Now run full workflow with runner
    runner = WorkflowRunner(config)
    success, error = runner.run()
    
    if success:
        print("✓ Test passed!")
    else:
        print(f"✗ Test failed: {error}")


def example_4_batch_testing():
    """Example 4: Run multiple pipelines in sequence."""
    print("=" * 70)
    print("Example 4: Batch testing multiple pipelines")
    print("=" * 70)
    
    pipelines = [
        "mt_coverage_merge",
        # Add more pipelines as they're configured
    ]
    
    results = {}
    
    for pipeline_name in pipelines:
        print(f"\nTesting: {pipeline_name}")
        print("-" * 70)
        
        try:
            config = get_pipeline_config(pipeline_name, repo_root=".")
            runner = WorkflowRunner(config)
            success, error = runner.run()
            results[pipeline_name] = ("PASS" if success else "FAIL", error)
        except Exception as e:
            results[pipeline_name] = ("ERROR", str(e))
    
    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    for pipeline_name, (status, error) in results.items():
        symbol = "✓" if status == "PASS" else "✗"
        print(f"{symbol} {pipeline_name}: {status}")
        if error:
            print(f"  Error: {error}")


def example_5_skip_options():
    """Example 5: Using skip options for faster iteration."""
    print("=" * 70)
    print("Example 5: Using skip options")
    print("=" * 70)
    
    config = get_pipeline_config("mt_coverage_merge", repo_root=".")
    
    # Skip Docker build if image hasn't changed
    config.skip_build = True
    
    # Keep old run directories for inspection
    config.skip_cleanup = True
    
    # Force reset of fake-GCS data (useful if wrappers changed)
    config.reset_fake_gcs = True
    
    runner = WorkflowRunner(config)
    
    print("Running with options:")
    print(f"  skip_build: {config.skip_build}")
    print(f"  skip_cleanup: {config.skip_cleanup}")
    print(f"  reset_fake_gcs: {config.reset_fake_gcs}")
    print()
    
    success, error = runner.run()
    
    if success:
        print("✓ Test passed!")
    else:
        print(f"✗ Test failed: {error}")


if __name__ == "__main__":
    # Uncomment the example you want to run:
    
    # example_1_predefined_pipeline()
    # example_2_custom_pipeline()
    # example_3_full_control()
    # example_4_batch_testing()
    # example_5_skip_options()
    
    print("Examples available:")
    print("  1. example_1_predefined_pipeline() - Run predefined config")
    print("  2. example_2_custom_pipeline() - Run with custom config")
    print("  3. example_3_full_control() - Full control with individual managers")
    print("  4. example_4_batch_testing() - Batch test multiple pipelines")
    print("  5. example_5_skip_options() - Using skip options")
    print()
    print("To run an example, uncomment in the main block and run this file.")
