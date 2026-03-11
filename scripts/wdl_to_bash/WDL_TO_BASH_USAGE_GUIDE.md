# WDL to Bash Extractor - Usage Guide

## Quick Reference

### Most Common Command
```bash
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py <pipeline.wdl>
```
Creates `extracted_tasks/` directory with all extracted bash scripts.

### Test Mode (Check for Bash Errors)
```bash
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py <pipeline.wdl> --test
```
Scans for bash syntax errors in all tasks **without** creating files.

---

## Step-by-Step Tutorial

### Step 0 (Optional): Test for Bash Errors First
Before extracting, you can scan the WDL file for bash syntax issues:

```bash
# Test mode: identify syntax errors without extraction
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py all_of_us/mitochondria/pipeline.wdl --test
```

This will output a report showing which tasks have bash syntax errors. If errors are found, they should be fixed in the WDL file before extraction.

Example output:
```
======================================================================
WDL BASH ERROR DETECTION REPORT
======================================================================

Total Tasks Found: 12
Valid Tasks: 9
Tasks with Bash Errors: 3

TASKS WITH BASH SYNTAX ERRORS:

❌ build_vcf_shard_mt
   Error: bash: line 78: unexpected EOF while looking for matching `)'
   Source: /path/to/pipeline.wdl
```

### Step 1: Extract Tasks
```bash
# Extract from mitochondria pipeline (run from workspace root)
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py all_of_us/mitochondria/mitochondria_pipeline.wdl \
  -o my_extracted_tasks

# Review what was extracted
ls -la my_extracted_tasks/
```

Expected output:
```
total 500
-rw-r--r-- Filter_inputs.json
-rwxr-xr-x Filter.sh
-rw-r--r-- M2_inputs.json
-rwxr-xr-x M2.sh
... (40+ more tasks)
```

### Step 2: Review a Generated Script
```bash
# Look at a specific task script
cat my_extracted_tasks/M2.sh

# Check syntax (bash -n = check without running)
bash -n my_extracted_tasks/M2.sh
```

### Step 3: Understand Variable Requirements
```bash
# Extract variables from the script
grep "TODO" my_extracted_tasks/M2.sh
```

Output shows required variables:
```bash
input_bam="${TODO_input_bam}"                # File - REQUIRED
ref_fasta="${TODO_ref_fasta}"               # File - REQUIRED
base_name="${TODO_base_name}"               # String - REQUIRED
```

### Step 4: Prepare Input Values
Option A - Use JSON template:
```bash
# Copy template
cp my_extracted_tasks/M2_inputs.json m2_inputs.json

# Edit with real values
nano m2_inputs.json
```

Edit to have:
```json
{
  "input_bam": "/data/samples/sample.bam",
  "ref_fasta": "/refs/mtDNA/rCRS.fasta",
  "base_name": "sample_001",
  ...
}
```

Option B - Edit script directly:
```bash
# Edit the script variables directly
nano my_extracted_tasks/M2.sh
```

Replace:
```bash
# Change this:
input_bam="${TODO_input_bam}"

# To this:
input_bam="/data/samples/sample.bam"
```

### Step 5: Execute the Script
```bash
# Method 1: Direct execution
bash my_extracted_tasks/M2.sh

# Method 2: Source and run
source my_extracted_tasks/M2.sh
M2_execution_function  # If no shebang, wrap in function

# Method 3: With explicit environment
export input_bam="/data/samples/sample.bam"
export ref_fasta="/refs/mtDNA/rCRS.fasta"
bash my_extracted_tasks/M2.sh
```

## Real-World Examples

### Example 1: Extract and Validate Quality Control Task
```bash
#!/bin/bash
set -euo pipefail

# Extract all tasks from pipeline
python3 wdl_to_bash_extractor.py pipelines/wdl/dna_seq/ExomeReprocessing.wdl \
  -o qc_scripts -l INFO

# Validate all bash scripts
echo "Validating extracted scripts..."
for script in qc_scripts/*.sh; do
  if bash -n "$script" 2>/dev/null; then
    echo "✓ $(basename $script)"
  else
    echo "✗ $(basename $script) - SYNTAX ERROR"
  fi
done

# Count extraction results
echo "Total tasks extracted: $(ls qc_scripts/*.sh | wc -l)"
```

### Example 2: Batch Extract Multiple Pipelines
```bash
#!/bin/bash
set -euo pipefail

# Process all pipelines in WARP
pipelines=(
  "pipelines/wdl/dna_seq/GermlineSingleSample.wdl"
  "pipelines/wdl/rna_seq/Optimus.wdl"
  "all_of_us/mitochondria/mitochondria_pipeline.wdl"
)

for pipeline in "${pipelines[@]}"; do
  pipeline_name=$(basename "${pipeline%.wdl}")
  echo "Extracting: $pipeline_name"
  
  python3 wdl_to_bash_extractor.py "$pipeline" \
    -o "extracted/$pipeline_name" \
    --no-input-templates  # Skip JSON for this batch
done

# Summary
echo ""
echo "=== EXTRACTION SUMMARY ==="
for dir in extracted/*/; do
  task_count=$(ls "$dir"*.sh 2>/dev/null | wc -l)
  printf "%-30s: %3d tasks\n" "$(basename $dir)" "$task_count"
done
```

### Example 3: Set Up Parameters for Docker Execution
```bash
#!/bin/bash
# Setup a task for Docker execution with mounted volumes

TASK_NAME="Filter"
SCRIPT_DIR="extracted_tasks"
WORK_DIR="/tmp/task_$$"
DATA_DIR="/data/analysis"

# Create working directory
mkdir -p "$WORK_DIR"

# Source the extracted script
source "$SCRIPT_DIR/${TASK_NAME}.sh"

# Override variables for docker context
raw_vcf="$DATA_DIR/calls/raw_variants.vcf"
ref_fasta="$DATA_DIR/reference/hg38.fasta"
ref_fai="$DATA_DIR/reference/hg38.fasta.fai"
ref_dict="$DATA_DIR/reference/hg38.dict"
base_name="sample_v1"
max_alt_allele_count=3
compress=false
run_contamination=true

# Export for task script
export raw_vcf ref_fasta ref_fai ref_dict base_name max_alt_allele_count compress run_contamination

# Run in docker with mounted volumes
docker run --rm \
  -v "$DATA_DIR:$DATA_DIR:ro" \
  -v "$WORK_DIR:$WORK_DIR" \
  -w "$WORK_DIR" \
  gatk:latest \
  bash "$SCRIPT_DIR/${TASK_NAME}.sh"
```

### Example 4: Create Pipeline DAG from Extracted Tasks
```bash
#!/bin/bash
# Generate a makefile from extracted tasks (for GNU make)

cat > Makefile.auto << 'EOF'
.PHONY: all clean validate

TASK_DIR := extracted_tasks
TASKS := $(basename $(notdir $(wildcard $(TASK_DIR)/*.sh)))

all: $(TASKS)

$(TASKS): 
	@echo "Running task: $@"
	bash $(TASK_DIR)/$@.sh

validate:
	@for script in $(TASK_DIR)/*.sh; do \
		if ! bash -n "$$script"; then \
			echo "Syntax error in $$script"; \
			exit 1; \
		fi; \
	done
	@echo "All scripts validated successfully"

clean:
	rm -rf $(TASK_DIR)/*.log *.sam *.bam

.DEFAULT_GOAL := all
EOF

echo "Generated Makefile.auto - run 'make' to execute tasks"
```

### Example 5: Monitor Execution with Logging
```bash
#!/bin/bash
# Execute extracted task with comprehensive logging

TASK_NAME="M2"
SCRIPT_PATH="extracted_tasks/${TASK_NAME}.sh"
LOG_FILE="logs/${TASK_NAME}_$(date +%Y%m%d_%H%M%S).log"

mkdir -p logs

# Run with logging
{
  echo "=== Execution Log for $TASK_NAME ==="
  echo "Start time: $(date)"
  echo "Script: $SCRIPT_PATH"
  echo "========================================"
  echo ""
  
  # Execute with verbose output
  bash -x "$SCRIPT_PATH" 2>&1
  
  exit_code=$?
  echo ""
  echo "========================================"
  echo "Exit code: $exit_code"
  echo "End time: $(date)"
} | tee "$LOG_FILE"

exit $exit_code
```

## Advanced Scenarios

### Scenario: Using Extracted Tasks in CI/CD Pipeline
```yaml
# GitHub Actions example
name: Run Extracted WDL Tasks

on: [push]

jobs:
  extract-and-test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      
      - name: Install miniwdl
        run: pip install miniwdl
      
      - name: Extract WDL tasks
        run: |
          python3 wdl_to_bash_extractor.py \
            pipelines/wdl/dna_seq/GermlineSingleSample.wdl \
            -o extracted_tasks
      
      - name: Validate extracted scripts
        run: |
          for script in extracted_tasks/*.sh; do
            bash -n "$script" || exit 1
          done
      
      - name: Run task tests
        run: |
          bash extracted_tasks/Align.sh
          bash extracted_tasks/MarkDuplicates.sh
```

### Scenario: Parallel Execution with GNU Parallel
```bash
#!/bin/bash
# Execute multiple tasks in parallel

# Extract all tasks
python3 wdl_to_bash_extractor.py pipeline.wdl -o tasks

# Run in parallel (4 at a time)
ls tasks/*.sh | \
  parallel -j 4 'echo "Running {}" && bash {}'

# With error reporting
ls tasks/*.sh | \
  parallel -j 4 --halt soon,fail=1 \
  'bash {} || echo "FAILED: {}"'
```

## Debugging Tips

### Debug Import Resolution
```bash
# See which imports were found and resolved
python3 wdl_to_bash_extractor.py pipeline.wdl -l DEBUG 2>&1 | grep -i import
```

### Debug Variable Substitution
```bash
# See all variable substitutions performed
python3 wdl_to_bash_extractor.py pipeline.wdl -l DEBUG 2>&1 | grep -i substit
```

### Verify Script Contents
```bash
# Count extracted tasks
ls extracted_tasks/*.sh | wc -l

# Show variable requirements
grep -h "TODO" extracted_tasks/*.sh | sort -u

# Find scripts with specific commands
grep -l "gatk" extracted_tasks/*.sh
```

### Test Individual Commands
```bash
# Extract just the command portion
tail -20 extracted_tasks/TaskName.sh

# Test command with dry-run (if command supports it)
bash -n extracted_tasks/TaskName.sh && echo "Syntax OK"
```

## Common Issues & Solutions

### Issue: "WDL module not found"
```bash
# Solution: Install miniwdl
pip install miniwdl --user
```

### Issue: "Could not resolve import"
```bash
# Check import paths in WDL file
grep "^import" pipelines/wdl/task.wdl

# Current directory must be repository root
cd /Users/juyang/codes/misc/warp
python3 wdl_to_bash_extractor.py ...
```

### Issue: "Bash syntax error" in generated script
```bash
# Review the command block for WDL-specific syntax
cat extracted_tasks/TaskName.sh | tail -30

# May need manual conversion of WDL expressions
# Example: ~{default="/path" variable} → ${variable:-/path}
```

### Issue: Variables showing as literal text
```bash
# Check if variables were properly substituted
grep -E '\$\{TODO_' extracted_tasks/*.sh | wc -l

# If count is high, check for script validation failures
python3 -c "
import subprocess
script = 'extracted_tasks/TaskName.sh'
result = subprocess.run(['bash', '-n'], stdin=open(script), capture_output=True)
print(result.stderr.decode())
"
```

## Performance Optimization

### For Large Pipelines
```bash
# Skip validation for faster extraction
python3 wdl_to_bash_extractor.py large_pipeline.wdl --no-validation

# Skip input JSON generation if not needed
python3 wdl_to_bash_extractor.py large_pipeline.wdl --no-input-templates
```

### For Batch Processing
```bash
# Create a job submission script
for wdl in pipelines/wdl/**/*.wdl; do
  dir=$(dirname "$wdl")
  name=$(basename "${wdl%.wdl}")
  
  # Queue extraction job (adjust for your job scheduler)
  sbatch -J "extract_${name}" <<EOF
#!/bin/bash
python3 wdl_to_bash_extractor.py "$wdl" -o "extracted/${name}"
EOF
done
```

## Next Steps

1. **Review extracted scripts** - Understand what commands will be executed
2. **Prepare input files** - Gather all required input data
3. **Test with small datasets** - Validate before full pipeline runs
4. **Integrate with automation** - Use in CI/CD, cron jobs, or job schedulers
5. **Monitor execution** - Track success/failures with logging

## Getting Help

For issues related to:
- **miniwdl parsing** → Check WDL syntax in original files
- **Variable conversion** → Review TODO markers in generated scripts  
- **Command execution** → Prepare input variables and test with `bash -n`
- **Import resolution** → Run with `-l DEBUG` flag

See `WDL_TO_BASH_README.md` for comprehensive documentation.
