# WARP Dev Container Setup - Summary

**Created on:** March 11, 2026

## Overview

A complete VS Code Dev Container configuration has been created for the WARP repository. This enables developers to work in a fully containerized environment with all necessary tools and dependencies pre-configured.

## Files Created

### Core Dev Container Configuration

**`.devcontainer/devcontainer.json`** (Primary Config)
- VS Code dev container definition
- Configures Python, Docker, and development tools
- Sets up useful aliases and helper commands
- Includes VS Code extensions for WDL/Python development
- Enables port forwarding for web applications

**`.devcontainer/Dockerfile`** (Container Image)
- Builds on Ubuntu 22.04 base image
- Installs Java 17, Python 3, Git, and build tools
- Pre-installs all Python dependencies
- Sets up womtool for WDL validation
- Optimized layers for faster builds

### Setup and Initialization Scripts

**`.devcontainer/onCreateCommand.sh`** (System Setup)
- Runs during container build
- Updates system packages
- Installs Java, Python, Git
- Downloads and configures womtool (WDL validator)
- ~2-3 minutes to complete

**`.devcontainer/postCreateCommand.sh`** (Python Dependencies)
- Runs after container creation
- Installs Python packages for WDL, testing, and development
- Sets up helpful bash aliases and functions
- Installs local requirements files if present
- ~1-2 minutes to complete

**`.devcontainer/postStartCommand.sh`** (Welcome Message)
- Runs each time the container starts
- Displays welcome message with quick-start commands
- Shows available documentation links

**`.devcontainer/setup-dev-container.sh`** (Interactive Setup)
- Optional manual setup script
- Validates environment and tools
- Downloads missing components
- Shows quick command reference

**`.devcontainer/download-womtool.sh`** (Utility)
- Helper script to download womtool JAR separately
- Useful if download fails during container build
- Supports different womtool versions

### Documentation

**`.devcontainer/README.md`** (Main Guide)
- Comprehensive setup and usage guide
- Quick start instructions (3 methods)
- List of all installed tools and packages
- Common commands with examples
- Troubleshooting section
- Performance tips and best practices

**`.devcontainer/codespaces.yml`** (GitHub Codespaces Config)
- Enables GitHub Codespaces support
- Prebuild configuration for faster startup
- Defines trigger paths for automated prebuild

### Docker Composition

**`docker-compose.yml`** (Container Orchestration)
- Enables easy docker-compose workflow
- Two services: warp-dev and jupyter
- Automatic volume mounting
- Port forwarding configuration
- Network setup for service communication

**`.devcontainer/.gitignore`** (Repository Filter)
- Prevents container artifacts from being committed
- Excludes womtool JAR files
- Excludes Python caches and build artifacts
- Excludes IDE and testing output

## Installed Tools & Packages

### System Tools
- **Java 17** - WDL validation with womtool
- **Python 3** - Latest Python version
- **Git** - Version control
- **Build Tools** - Development utilities

### Python Packages (24 packages)

**WDL & Testing:**
- miniwdl 1.9.1 - WDL execution/validation
- pytest 7.4.3 - Testing framework
- pytest-xdist 3.5.0 - Parallel testing
- pytest-timeout 2.2.0 - Test timeouts
- pytest-cov 4.1.0 - Coverage reporting

**Code Quality:**
- black 23.12.0 - Code formatter
- ruff 0.1.11 - Fast Python linter
- pylint 3.0.3 - Code analysis
- mypy 1.8.0 - Static type checking

**Data & Templates:**
- pyyaml 6.0.1 - YAML processing
- jinja2 3.1.6 - Template engine
- jsonschema 4.20.0 - JSON validation

**Interactive Development:**
- IPython - Enhanced Python shell
- Jupyter - Interactive notebooks
- JupyterLab - Modern notebook interface

### VS Code Extensions (9 extensions)
- Python language support
- Pylance type checking
- Ruff linter integration
- Black formatter
- Docker support
- GitHub Copilot
- YAML support
- And more...

## Bash Aliases & Functions

The following aliases are automatically available inside the container:

```bash
# WDL Validation
validate-wdl <file>              # Quick validation
validate-all-wdl [dir]           # Validate directory

# Extraction
extract-bash <pipeline.wdl>      # Extract to bash scripts

# Testing
run-tests                        # Run pytest suite

# Information
warp-info                        # Repository overview
warp-version                     # Show pipeline versions
```

## Getting Started

### Method 1: VS Code (Recommended)
1. Open WARP repository in VS Code
2. Install "Dev Containers" extension
3. Click >< button in bottom-left corner
4. Select "Reopen in Container"
5. Wait for build (~5-10 minutes first time)

### Method 2: Docker Compose
```bash
docker-compose up -d
docker-compose exec warp-dev bash
```

### Method 3: Manual Docker
```bash
docker build -f .devcontainer/Dockerfile -t warp-dev:latest .
docker run -it -v $(pwd):/workspace warp-dev:latest bash
```

## Key Features

✅ **All-in-One Setup** - No manual tool installation needed
✅ **Automatic Dependency Management** - Python and system packages pre-installed
✅ **VS Code Integration** - Seamless development experience with full extension support
✅ **Port Forwarding** - Access Jupyter and applications on localhost
✅ **Git Configuration** - Automatic git config and SSH key mounting
✅ **Docker Support** - Docker-in-Docker enabled for container builds
✅ **Performance Optimized** - Cached volumes for better file access
✅ **Reproducible Environment** - Same setup across all developers

## Common Workflows

### Validate WDL Files
```bash
womtool validate pipelines/wdl/some_pipeline/SomePipeline.wdl
```

### Run Tests with Coverage
```bash
python3 -m pytest --cov -v all_of_us/mitochondria/test_local/
```

### Extract Bash Scripts
```bash
python3 wdl_to_bash_extractor.py pipelines/wdl/some_pipeline/SomePipeline.wdl -o ./extracted_tasks
```

### Use Jupyter Lab
```bash
jupyter lab --ip=0.0.0.0 --allow-root
# Access at http://localhost:8888
```

## Troubleshooting

If womtool not found:
```bash
bash .devcontainer/download-womtool.sh
```

If packages missing:
```bash
python3 -m pip install --upgrade --force-reinstall miniwdl pytest pytest-xdist ...
```

To rebuild container:
- VS Code: Ctrl+Shift+P → "Dev Containers: Rebuild Container"
- Docker Compose: `docker-compose down -v && docker-compose up --build -d`

## Project Structure

```
.devcontainer/
├── devcontainer.json           # Main VS Code config
├── Dockerfile                  # Container image definition  
├── onCreateCommand.sh          # System setup (runs during build)
├── postCreateCommand.sh        # Python setup (runs after create)
├── postStartCommand.sh         # Welcome message (runs on start)
├── setup-dev-container.sh      # Interactive setup utility
├── download-womtool.sh         # womtool download helper
├── codespaces.yml              # GitHub Codespaces config
├── .gitignore                  # Container artifacts filter
└── README.md                   # Complete user guide

docker-compose.yml              # Optional: Docker Compose orchestration
```

## Next Steps

1. **Review** `.devcontainer/README.md` for detailed documentation
2. **Open Repository** in VS Code and reopen in container
3. **Run** `show-versions` to verify everything is installed
4. **Validate** a WDL file: `womtool validate <file.wdl>`
5. **Start Developing!**

## Notes

- Container runs as root (appropriate for development)
- All path changes are immediately reflected on host
- Python environment is preserved across sessions
- Extension settings sync with VS Code account if configured
- GitHub Codespaces support enabled for remote development

## Documentation Links

- WDL Specification: https://openwdl.org/
- miniwdl GitHub: https://github.com/chanzuckerberg/miniwdl
- Dev Containers: https://containers.dev/
- Docker Documentation: https://docs.docker.com/

---

**Setup completed successfully!** You can now develop on WARP with a fully containerized environment.
