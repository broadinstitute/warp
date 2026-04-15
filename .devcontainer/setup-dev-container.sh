#!/bin/bash
# setup-dev-container.sh - Quick setup script for WARP dev container
# This script helps set up the dev container environment quickly

set -e

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}═════════════════════════════════════════${NC}"
echo -e "${BLUE}  WARP Dev Container Setup Script${NC}"
echo -e "${BLUE}═════════════════════════════════════════${NC}"
echo ""

# Check if running in dev container
if [ ! -f "/.dockerenv" ]; then
  echo -e "${YELLOW}⚠️  This script should be run inside the dev container${NC}"
  echo ""
  echo "To enter the dev container:"
  echo "  1. Press Ctrl+Shift+P in VS Code"
  echo "  2. Type 'Dev Containers: Reopen in Container'"
  echo ""
  echo "Or use Docker Compose:"
  echo "  docker-compose up -d"
  echo "  docker-compose exec warp-dev /bin/bash"
  echo ""
  exit 0
fi

echo -e "${GREEN}✓ Running inside dev container${NC}"
echo ""

# Check environment
echo -e "${BLUE}Checking environment...${NC}"
echo "  Java:   $(java -version 2>&1 | tail -n2 | head -n1)"
echo "  Python: $(python3 --version)"
echo "  Git:    $(git --version)"
echo ""

# Check womtool
echo -e "${BLUE}Checking womtool...${NC}"
if command -v womtool &> /dev/null; then
  echo -e "  ${GREEN}✓ womtool is available${NC}"
  womtool --version 2>/dev/null || echo "  Version check completed"
else
  echo -e "  ${YELLOW}⚠️  womtool not found${NC}"
  echo "  Downloading womtool..."
  bash .devcontainer/download-womtool.sh
fi
echo ""

# Check Python packages
echo -e "${BLUE}Checking Python packages...${NC}"
python3 -c "
import miniwdl
import pytest
import yaml
import jinja2
print('  ✓ Core packages installed')
" || {
  echo -e "  ${YELLOW}⚠️  Some packages missing, installing...${NC}"
  python3 -m pip install --quiet miniwdl pytest pyyaml jinja2 jsonschema
}
echo ""

# Display quick commands
echo -e "${GREEN}═════════════════════════════════════════${NC}"
echo -e "${GREEN}Setup Complete! Ready to develop.${NC}"
echo -e "${GREEN}═════════════════════════════════════════${NC}"
echo ""

echo "Quick Start Commands:"
echo ""
echo -e "${BLUE}Validation:${NC}"
echo "  womtool validate <file.wdl>         # Validate WDL syntax"
echo "  validate-all-wdl <directory>        # Validate all WDL files"
echo ""

echo -e "${BLUE}Extraction:${NC}"
echo "  python3 wdl_to_bash_extractor.py <pipeline.wdl>"
echo ""

echo -e "${BLUE}Testing:${NC}"
echo "  python3 -m pytest -v                # Run all tests"
echo "  python3 -m pytest -n auto          # Run tests in parallel"
echo ""

echo -e "${BLUE}Information:${NC}"
echo "  warp-info                           # Show repo overview"
echo "  show-versions                       # Show tool versions"
echo ""

echo "For more details, see: .devcontainer/README.md"
echo ""
