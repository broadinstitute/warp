#!/bin/bash
# postCreateCommand.sh - Runs after container is created to install Python dependencies
set -e

echo "📚 Installing Python dependencies..."

# Install core Python packages for WDL development
python3 -m pip install --quiet \
  miniwdl==1.9.1 \
  pytest==7.4.3 \
  pytest-xdist==3.5.0 \
  pytest-timeout==2.2.0 \
  pytest-cov==4.1.0 \
  pyyaml==6.0.1 \
  jinja2==3.1.6 \
  jsonschema==4.20.0 \
  black==23.12.0 \
  pylint==3.0.3 \
  mypy==1.8.0 \
  ruff==0.1.11

# Install development tools
python3 -m pip install --quiet \
  ipython \
  jupyter \
  jupyterlab

# Install any local requirements files if they exist
echo "🔍 Checking for local requirements files..."

if [ -f "requirements.txt" ]; then
  echo "Installing root requirements.txt..."
  python3 -m pip install --quiet -r requirements.txt
fi

if [ -f "all_of_us/mitochondria/requirements.txt" ]; then
  echo "Installing mitochondria test requirements..."
  python3 -m pip install --quiet -r all_of_us/mitochondria/requirements.txt 2>/dev/null || true
fi

if [ -f "verification/test-wdls/scripts/requirements.txt" ]; then
  echo "Installing verification test requirements..."
  python3 -m pip install --quiet -r verification/test-wdls/scripts/requirements.txt 2>/dev/null || true
fi

# Create useful aliases and helper commands
cat >> /root/.bashrc << 'EOF'

# WARP Development Helpers
alias validate-wdl='womtool validate'
alias extract-bash='python3 wdl_to_bash_extractor.py'
alias run-tests='python3 -m pytest -v'
alias warp-version='cat pipeline_versions.txt | head -20'

# Function to validate all WDL files in a directory
validate-all-wdl() {
  local dir=${1:-.}
  echo "Validating all WDL files in $dir..."
  find "$dir" -name "*.wdl" -type f | while read -r wdl; do
    echo "Validating: $wdl"
    womtool validate "$wdl" || echo "❌ Failed: $wdl"
  done
}

# Function to show WARP directory overview
warp-info() {
  echo "🧬 WARP Repository Overview"
  echo "============================"
  echo ""
  echo "Pipelines:"
  ls -d pipelines/wdl/*/ 2>/dev/null | xargs -I {} basename {} | head -10
  echo ""
  echo "Tasks:"
  ls tasks/wdl/*.wdl 2>/dev/null | wc -l | xargs echo "Total task files:"
  echo ""
  echo "Structures:"
  find structs -name "*.wdl" 2>/dev/null | wc -l | xargs echo "Total struct files:"
  echo ""
  echo "Latest Pipeline Versions:"
  head -5 pipeline_versions.txt 2>/dev/null || echo "pipeline_versions.txt not found"
}
EOF

echo "✅ Python dependencies installed successfully!"
echo ""
echo "📋 Available commands:"
echo "  - womtool validate <file.wdl>   : Validate WDL files"
echo "  - validate-all-wdl [dir]        : Validate all WDL files in directory"
echo "  - validate-wdl <file.wdl>       : Alias for womtool validate"
echo "  - extract-bash <pipeline.wdl>   : Extract bash scripts from WDL"
echo "  - run-tests                     : Run pytest suite"
echo "  - warp-info                     : Show WARP repository overview"
echo "  - warp-version                  : Show pipeline versions"
echo ""
