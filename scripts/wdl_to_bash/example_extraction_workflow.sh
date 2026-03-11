#!/bin/bash
# Example: Complete workflow for extracting and executing WDL tasks
# This demonstrates the full pipeline from extraction through execution

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACTION_DIR="${SCRIPT_DIR}/extracted_tasks_example"
WORK_DIR="${SCRIPT_DIR}/task_work_$$"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo_info() {
  echo -e "${BLUE}[INFO]${NC} $1"
}

echo_success() {
  echo -e "${GREEN}[SUCCESS]${NC} $1"
}

echo_warn() {
  echo -e "${YELLOW}[WARN]${NC} $1"
}

echo_error() {
  echo -e "${RED}[ERROR]${NC} $1"
}

cleanup() {
  echo_info "Cleaning up temporary directory: $WORK_DIR"
  rm -rf "$WORK_DIR"
}

trap cleanup EXIT

main() {
  echo_info "Starting WDL Task Extraction and Execution Example"
  echo_info "===================================================="
  echo ""
  
  # Step 1: Extract tasks from mitochondria pipeline
  echo_info "Step 1: Extracting WDL tasks..."
  
  # Determine workspace root based on script location
  WORKSPACE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
  
  if [ ! -f "${WORKSPACE_ROOT}/all_of_us/mitochondria/mitochondria_pipeline.wdl" ]; then
    echo_error "Cannot find mitochondria_pipeline.wdl"
    echo_error "Expected at: ${WORKSPACE_ROOT}/all_of_us/mitochondria/mitochondria_pipeline.wdl"
    return 1
  fi
  
  python3 "${SCRIPT_DIR}/wdl_to_bash_extractor.py" \
    "${WORKSPACE_ROOT}/all_of_us/mitochondria/mitochondria_pipeline.wdl" \
    -o "$EXTRACTION_DIR" \
    -l INFO
  
  echo_success "Task extraction completed"
  echo ""
  
  # Step 2: Validate all extracted scripts
  echo_info "Step 2: Validating extracted bash scripts..."
  
  total_scripts=$(ls "$EXTRACTION_DIR"/*.sh 2>/dev/null | wc -l)
  valid_scripts=0
  invalid_scripts=()
  
  for script in "$EXTRACTION_DIR"/*.sh; do
    if bash -n "$script" 2>/dev/null; then
      ((valid_scripts++))
    else
      invalid_scripts+=("$(basename $script)")
    fi
  done
  
  echo_success "Validation complete: $valid_scripts/$total_scripts scripts are valid"
  
  if [ ${#invalid_scripts[@]} -gt 0 ]; then
    echo_warn "Invalid scripts found:"
    for script in "${invalid_scripts[@]}"; do
      echo "  - $script"
    done
  fi
  echo ""
  
  # Step 3: Show sample task information
  echo_info "Step 3: Analyzing extracted tasks..."
  
  # Count tasks by category
  echo_success "Generated task files:"
  find "$EXTRACTION_DIR" -name "*.sh" | head -10 | xargs basename -a
  echo "  ... and $(ls $EXTRACTION_DIR/*.sh | wc -l) total"
  echo ""
  
  # Show input requirements for M2 task
  if [ -f "$EXTRACTION_DIR/M2.sh" ]; then
    echo_info "Step 3a: Example - M2 Task Requirements"
    echo "Input variables for M2 task:"
    grep "TODO" "$EXTRACTION_DIR/M2.sh" | sed 's/^/  /'
  fi
  echo ""
  
  # Step 4: Show task inputs template
  echo_info "Step 4: Task Input Templates"
  
  if [ -f "$EXTRACTION_DIR/Filter_inputs.json" ]; then
    echo_success "Filter task input template:"
    echo "$(cat $EXTRACTION_DIR/Filter_inputs.json | head -10)"
    echo "  ..."
  fi
  echo ""
  
  # Step 5: Demonstrate task execution with placeholder values
  echo_info "Step 5: Demonstrating Task Execution Pattern"
  
  # Create work directory
  mkdir -p "$WORK_DIR"
  
  # Copy a simple task script
  if [ -f "$EXTRACTION_DIR/ChainSwap.sh" ]; then
    echo_success "Using ChainSwap task as example"
    
    # Show task contents
    echo "Task source code:"
    echo "---"
    head -30 "$EXTRACTION_DIR/ChainSwap.sh" | sed 's/^/  /'
    echo "  ..."
    echo "---"
    echo ""
    
    # Prepare execution environment
    echo_info "Setting up execution environment for ChainSwap..."
    
    # Create a wrapper script that demonstrates execution
    cat > "$WORK_DIR/demo_execution.sh" << 'EOF'
#!/bin/bash
# Demo script showing how to execute extracted tasks

# Define input variables (these would come from actual files in real execution)
chain_file="/path/to/chain.txt"
query_file="/path/to/query.bed"
target_file="/path/to/target.bed" 
output_prefix="demo_output"

echo "Would execute ChainSwap with:"
echo "  chain_file: $chain_file"
echo "  query_file: $query_file"
echo "  target_file: $target_file"
echo "  output_prefix: $output_prefix"

# In real execution, uncomment:
# bash "$EXTRACTION_DIR/ChainSwap.sh"
EOF
    
    chmod +x "$WORK_DIR/demo_execution.sh"
    
    echo_success "Execution wrapper created at: $WORK_DIR/demo_execution.sh"
  fi
  echo ""
  
  # Step 6: Generate summary report
  echo_info "Step 6: Extraction Summary"
  echo "============================"
  echo ""
  echo "Location: $EXTRACTION_DIR"
  echo "Total scripts: $(ls $EXTRACTION_DIR/*.sh 2>/dev/null | wc -l)"
  echo "Total input templates: $(ls $EXTRACTION_DIR/*_inputs.json 2>/dev/null | wc -l)"
  echo ""
  
  # Count by task type (heuristic)
  echo "Sample tasks extracted:"
  ls "$EXTRACTION_DIR"/*.sh | head -5 | xargs -I {} basename {} .sh | sed 's/^/  - /'
  
  local remaining=$(($(ls $EXTRACTION_DIR/*.sh | wc -l) - 5))
  if [ $remaining -gt 0 ]; then
    echo "  ... and $remaining more tasks"
  fi
  echo ""
  
  # Step 7: Usage instructions
  echo_info "Next Steps"
  echo "=========="
  echo ""
  echo "1. Review generated scripts:"
  echo "   $ cat $EXTRACTION_DIR/TaskName.sh"
  echo ""
  echo "2. Check input requirements:"
  echo "   $ grep 'TODO' $EXTRACTION_DIR/*.sh"
  echo ""
  echo "3. Create input files:"
  echo "   $ cp $EXTRACTION_DIR/TaskName_inputs.json task_inputs.json"
  echo "   $ nano task_inputs.json  # Edit with real values"
  echo ""
  echo "4. Execute a task:"
  echo "   $ bash $EXTRACTION_DIR/TaskName.sh"
  echo ""
  echo "5. Monitor execution:"
  echo "   $ bash -x $EXTRACTION_DIR/TaskName.sh 2>&1 | tee execution.log"
  echo ""
  
  echo_success "Example workflow completed!"
  echo ""
  echo "For more information, see: WDL_TO_BASH_README.md"
}

# Run main function
main "$@"
