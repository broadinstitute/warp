# WARP Dev Container Setup

This directory contains the VS Code Dev Container configuration for the WARP (Warp Analysis Research Pipelines) repository.

## Quick Start

### Option 1: VS Code with Remote Containers Extension (Recommended)

1. **Install Extensions**
   - Install the "Dev Containers" extension in VS Code
   - Install the "Remote - Containers" extension

2. **Open in Container**
   - Open the WARP repository in VS Code
   - Click the green button in the bottom-left corner (><)
   - Select "Reopen in Container"
   - Wait for the container to build and start

3. **Done!**
   - You're now in the development environment
   - All tools and dependencies are pre-installed
   - Terminal commands are ready to use

### Option 2: Docker Compose

```bash
# Build and start the development container
docker-compose up -d

# Enter the container shell
docker-compose exec warp-dev /bin/bash

# For Jupyter Lab
docker-compose exec jupyter bash
# Access at http://localhost:8889
```

### Option 3: Manual Docker Build

```bash
# Build the container
docker build -f .devcontainer/Dockerfile -t warp-dev:latest .

# Run the container
docker run -it -v $(pwd):/workspace warp-dev:latest bash
```

## What's Installed

### System Tools
- **Java 17** - Required for womtool WDL validation
- **Python 3** - Latest Python 3 version
- **Git** - Version control
- **Build Tools** - Essential development utilities
- **Docker** - For containerized workflows (Docker-in-Docker supported)

### Python Packages

#### WDL Development
- **miniwdl** (1.9.1) - WDL execution and validation
- **pytest** (7.4.3) - Testing framework
- **pytest-xdist** - Parallel test execution
- **pytest-timeout** - Test timeouts
- **pytest-cov** - Code coverage reporting

#### Code Quality
- **black** (23.12.0) - Code formatter
- **ruff** (0.1.11) - Fast Python linter
- **pylint** (3.0.3) - Code analysis
- **mypy** (1.8.0) - Static type checker

#### Data & Templates
- **pyyaml** (6.0.1) - YAML processing
- **jinja2** (3.1.6) - Template engine
- **jsonschema** (4.20.0) - JSON validation

#### Interactive Development
- **IPython** - Enhanced Python shell
- **Jupyter** - Interactive notebooks
- **JupyterLab** - Modern Jupyter interface

### VS Code Extensions
- **Python** - Python language support
- **Pylance** - Python type checking
- **Ruff** - Linter integration
- **Black Formatter** - Code formatting
- **Docker** - Docker support
- **GitHub Copilot** - AI assistance
- **YAML** - YAML file support
- **Git** - Git integration

## Common Commands

### WDL Validation

```bash
# Validate a single WDL file
womtool validate pipelines/wdl/some_pipeline/SomePipeline.wdl

# Validate all WDL files in a directory
validate-all-wdl pipelines/wdl/mitochondria/

# Quick alias
validate-wdl pipeline.wdl
```

### WDL Extraction

```bash
# Extract bash scripts from WDL pipeline
python3 wdl_to_bash_extractor.py pipelines/wdl/some_pipeline/SomePipeline.wdl

# Extract with custom output directory
python3 wdl_to_bash_extractor.py pipeline.wdl -o ./my_tasks

# Extract with debug logging
python3 wdl_to_bash_extractor.py pipeline.wdl -l DEBUG 2>&1 | tee debug.log

# Quick alias
extract-bash pipeline.wdl
```

### Testing

```bash
# Run all tests
python3 -m pytest -v

# Run tests in a specific directory
python3 -m pytest -v all_of_us/mitochondria/test_local/

# Run tests with coverage
python3 -m pytest --cov=. -v

# Run tests in parallel
python3 -m pytest -n auto -v

# Quick alias
run-tests
```

### Repository Info

```bash
# Show quick repository overview
warp-info

# Display pipeline versions
warp-version

# Check environment versions
show-versions
```

## File Structure

```
.devcontainer/
├── devcontainer.json       # Main VS Code dev container config
├── Dockerfile              # Docker image definition
├── onCreateCommand.sh      # Runs during container build (system setup)
├── postCreateCommand.sh    # Runs after container creation (Python deps)
├── postStartCommand.sh     # Runs each container start (welcome message)
├── download-womtool.sh     # Helper to download womtool JAR
└── README.md              # This file

docker-compose.yml          # Docker Compose configuration (optional)
```

## Environment Variables

Inside the container:

```bash
WARP_HOME=/workspace              # Root of the repository
DEBIAN_FRONTEND=noninteractive    # Non-interactive package installs
PYTHONUNBUFFERED=1               # Python output buffering disabled
```

## Configuring Python Path

The dev container automatically configures the Python interpreter. To verify:

```bash
which python3
# Output: /usr/bin/python3

python3 --version
# Output: Python 3.x.x

python3 -m pip list | head -10
```

## Troubleshooting

### womtool validation fails

If you get a "womtool not found" error:

```bash
# Download womtool manually
bash .devcontainer/download-womtool.sh

# Or specify the version
bash .devcontainer/download-womtool.sh 90
```

### Port conflicts

If ports 8080 or 8888 are already in use:

1. In `devcontainer.json`, change the `forwardPorts` values
2. Or in `docker-compose.yml`, modify the port mappings
3. Then rebuild the container

### Package installation issues

If Python packages fail to install:

```bash
# Reinstall all dependencies
python3 -m pip install --upgrade --force-reinstall \
  miniwdl pytest pytest-xdist pytest-timeout pytest-cov \
  pyyaml jinja2 jsonschema black pylint mypy ruff

# Clear pip cache
python3 -m pip cache purge
```

### Rebuild the container

If you need a clean rebuild:

```bash
# In VS Code: Press Ctrl+Shift+P and run "Dev Containers: Rebuild Container"

# Or via docker
docker-compose down -v
docker-compose build --no-cache
docker-compose up -d
```

## Git Configuration

Your git configuration is automatically mounted into the container:

```bash
git config --global user.name    # Shows your configured name
git config --global user.email   # Shows your configured email
```

SSH keys are also mounted if available in `~/.ssh/`.

## Port Forwarding

The dev container forwards these ports:

- **8080** - Application server port
- **8888** - Jupyter Lab port

Simply access them on your host machine:
- Application: http://localhost:8080
- Jupyter: http://localhost:8888

## Additional Documentation

Inside the repository:

- **README.md** - Main WARP documentation
- **WARP_WDL_Style_Guide.md** - WDL coding standards
- **WDL_TO_BASH_README.md** - WDL extraction tool details
- **QUICK_REFERENCE.md** - Quick reference card
- **.github/copilot-instructions.md** - Coding guidelines

## Performance Tips

1. **Use cached volumes** - Mounted volumes use `consistency=cached` for better performance
2. **Parallel testing** - Use `pytest -n auto` for faster test execution
3. **Quick validation** - Use `womtool validate` before commit
4. **IDE caching** - Python language server caches are preserved across sessions

## Support and Resources

- [Dev Containers Documentation](https://containers.dev/)
- [VS Code Remote Development](https://code.visualstudio.com/docs/remote/remote-overview)
- [WARP Repository](https://github.com/broadinstitute/warp)
- [WDL Specification](https://openwdl.org/)
- [miniwdl Documentation](https://github.com/chanzuckerberg/miniwdl)

## Notes

- The container runs as root, which is appropriate for development
- Docker-out-of-Docker is enabled if you need to build/run containers
- All changes to files are immediately reflected on your host machine
- The container persists Python packages and dependencies between sessions
