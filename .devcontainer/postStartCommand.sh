#!/bin/bash
# postStartCommand.sh - Runs each time the container starts

echo "🚀 Starting WARP development environment..."

# Display welcome message
echo ""
echo "🧬 Welcome to WARP Development Environment"
echo "=========================================="
echo ""
echo "📂 Repository: Warp Analysis Research Pipelines"
echo "🔗 Location: $WARP_HOME"
echo ""
echo "Quick Start Commands:"
echo "  • womtool validate <file.wdl>   - Validate WDL syntax"
echo "  • python3 wdl_to_bash_extractor.py <pipeline.wdl> - Extract tasks to bash"
echo "  • python3 -m pytest -v          - Run tests"
echo "  • warp-info                     - Show repository overview"
echo ""
echo "Documentation:"
echo "  • README.md                     - Main documentation"
echo "  • WDL_TO_BASH_README.md         - WDL extraction tool guide"
echo "  • WARP_WDL_Style_Guide.md       - WDL coding standards"
echo "  • QUICK_REFERENCE.md            - Quick reference card"
echo ""
