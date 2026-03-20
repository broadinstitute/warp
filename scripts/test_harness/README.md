# WARP Test Harness

A generalizable Python-based test framework for WARP WDL pipelines. Handles Docker image management, fake-GCS setup, WDL patching, and miniwdl workflow execution.

## Overview

The test harness automates the complete workflow testing process:

1. **Docker Management** - Build and shadow test Docker images
2. **GCS Setup** - Validate and populate fake-GCS server with test data
3. **WDL Patching** - Replace production Docker images with test versions
4. **Workflow Execution** - Run miniwdl with proper environment configuration
5. **Error Reporting** - Collect and display workflow errors

## Quick Start

### Requirements

- Docker and Docker Desktop (for fake-GCS server)
- miniwdl (for workflow execution)
- Python 3.7+
- fake-gcs-server running locally

### Start fake-GCS server

```bash
docker run -d -p 4443:4443 --name fake-gcs \
  fsouza/fake-gcs-server -scheme http \
  -external-url http://host.docker.internal:4443 \
  -public-host host.docker.internal:4443
```

### Run tests

```bash
# Using Python module
python -m scripts.test_harness mt_coverage_merge

# Using entry point script
./scripts/run_test.py mt_coverage_merge

# With options
./scripts/run_test.py mt_coverage_merge --skip-cleanup --reset-gcs
```

## Usage

### Predefined Pipelines

Available pipelines are configured in `config.py`:

```bash
# List available pipelines
./scripts/run_test.py --list

# Run specific pipeline
./scripts/run_test.py mt_coverage_merge
```

### Custom Pipelines

For any WDL pipeline not in the predefined list:

```bash
./scripts/run_test.py custom \
  --workflow pipelines/wdl/my_pipeline/pipeline.wdl \
  --inputs testing/my_inputs.json \
  --prod-image "us.gcr.io/my-org/image:v1.0" \
  --test-image "image:local-test" \
  --dockerfile testing/image/
```

### Command-line Options

```
--repo-root DIR          Repository root (default: current directory)
--skip-cleanup           Keep old run directories for inspection
--skip-build             Skip Docker image build (use existing image)
--reset-gcs              Force re-setup of fake-GCS data
--list                   Show available pipelines
--version                Show version
```

### Environment Variables

Alternative to command-line flags:

```bash
# Keep old run directories
WARP_SKIP_CLEANUP=true ./scripts/run_test.py mt_coverage_merge

# Skip Docker build
WARP_SKIP_BUILD=true ./scripts/run_test.py mt_coverage_merge

# Reset fake-GCS data
WARP_RESET_GCS=true ./scripts/run_test.py mt_coverage_merge

# Combine
WARP_SKIP_CLEANUP=true WARP_SKIP_BUILD=true ./scripts/run_test.py mt_coverage_merge
```

## Architecture

### Module Structure

```
scripts/test_harness/
├── __init__.py           # Package definition
├── __main__.py           # Module entry point
├── config.py             # Configuration management
├── docker_manager.py     # Docker operations
├── gcs_manager.py        # Fake-GCS validation and setup
├── wdl_patcher.py        # WDL file modifications
├── runner.py             # Workflow orchestration
└── cli.py                # Command-line interface
```

### Adding New Pipelines

1. Add configuration to `PIPELINE_CONFIGS` in `config.py`:

```python
PIPELINE_CONFIGS = {
    "my_pipeline": lambda repo_root: TestConfig(
        workflow_wdl="pipelines/wdl/my_pipeline/pipeline.wdl",
        inputs_json="testing/my_pipeline_inputs.json",
        docker=DockerConfig(
            prod_image="us.gcr.io/my-org/my-image:1.0.0",
            test_image="my-image:local-test",
            dockerfile_dir="testing/my-image/",
        ),
        repo_root=repo_root,
    ),
}
```

2. Run tests:

```bash
./scripts/run_test.py my_pipeline
```

### Key Classes

#### `TestConfig` (config.py)

Holds all test configuration including paths, Docker settings, and control flags.

**Usage:**
```python
from test_harness.config import TestConfig, DockerConfig

config = TestConfig(
    workflow_wdl="pipelines/wdl/my_pipeline/pipeline.wdl",
    inputs_json="testing/inputs.json",
    docker=DockerConfig(
        prod_image="us.gcr.io/my-org/image:v1.0",
        test_image="image:local-test",
        dockerfile_dir="testing/image/",
    ),
)
```

#### `DockerManager` (docker_manager.py)

Manages Docker image lifecycle: build, shadow (replace), restore.

**Usage:**
```python
from test_harness.docker_manager import DockerManager
from test_harness.config import DockerConfig

docker_config = DockerConfig(
    prod_image="us.gcr.io/my-org/image:v1.0",
    test_image="image:local-test",
    dockerfile_dir="/path/to/dockerfile",
)
manager = DockerManager(docker_config)

# Build test image
manager.build_test_image()

# Replace production image with test image
manager.shadow_prod_image()

# Restore (remove shadowing)
manager.restore_prod_image()
```

#### `GCSManager` (gcs_manager.py)

Manages fake-GCS server setup and test data validation.

**Usage:**
```python
from test_harness.gcs_manager import GCSManager
from test_harness.config import GCSConfig
from pathlib import Path

gcs_config = GCSConfig()
manager = GCSManager(gcs_config, Path("/repo/root"))

# Check if server is available
if manager.is_available():
    # Ensure test data exists, uploading if needed
    manager.ensure_test_data(force_reset=False)
```

#### `WDLPatcher` (wdl_patcher.py)

Patches WDL files to replace docker images for testing.

**Usage:**
```python
from test_harness.wdl_patcher import WDLPatcher
from pathlib import Path

patcher = WDLPatcher(Path("/repo/root"))

# Patch a WDL file
patcher.patch_workflow(
    workflow_path=Path("/repo/root/pipelines/wdl/my_pipeline/pipeline.wdl"),
    output_path=Path("/repo/root/tmp/pipeline.test.wdl"),
    image_replacements={
        "us.gcr.io/my-org/image:1.0.0": "image:local-test"
    }
)
```

#### `WorkflowRunner` (runner.py)

Orchestrates the complete test execution workflow.

**Usage:**
```python
from test_harness.runner import WorkflowRunner
from test_harness.config import get_pipeline_config

config = get_pipeline_config("mt_coverage_merge", "/repo/root")
runner = WorkflowRunner(config)

success, error = runner.run()
if success:
    print("Test passed!")
else:
    print(f"Test failed: {error}")
```

## Development

### Running as a module

```bash
# From repository root
python -m scripts.test_harness mt_coverage_merge
python -m scripts.test_harness --list
```

### Direct Python usage

```python
from scripts.test_harness.config import get_pipeline_config
from scripts.test_harness.runner import WorkflowRunner

config = get_pipeline_config("mt_coverage_merge", ".")
runner = WorkflowRunner(config)
success, error = runner.run()
```

## Troubleshooting

### fake-gcs-server not running

```
✗ fake-gcs-server is not running!

Start it with:
  docker run -d -p 4443:4443 --name fake-gcs \
    fsouza/fake-gcs-server -scheme http \
    -external-url http://host.docker.internal:4443 \
    -public-host host.docker.internal:4443
```

### miniwdl not found

Install miniwdl:
```bash
pip install miniwdl
```

### Docker permission denied

Add current user to docker group:
```bash
sudo usermod -aG docker $USER
```

### WDL validation errors

The test framework doesn't validate WDL before patching. Validate manually:
```bash
java -jar ~/womtool-90.jar validate pipelines/wdl/my_pipeline/pipeline.wdl
```

## Migration from Bash Script

The original `run_test.sh` is still available but deprecated. To migrate:

1. **Bash script** (old):
   ```bash
   ./run_test.sh
   ```

2. **Python entry point** (new):
   ```bash
   ./scripts/run_test.py mt_coverage_merge
   ```

3. **Python module** (alternative):
   ```bash
   python -m scripts.test_harness mt_coverage_merge
   ```

All three approaches produce identical results, but the Python-based test harness is more maintainable and extensible.
