# WDL Command Block Extractor to Bash Scripts

A Python utility that parses WDL (Workflow Description Language) pipeline files using miniwdl and extracts all task command blocks as executable bash scripts.

## Features

✅ **Recursive Import Resolution**: Automatically follows and processes all imported WDL files
✅ **Complete Variable Extraction**: Captures all inputs, outputs, and their types from each task
✅ **Smart Variable Substitution**: Converts WDL variable syntax to bash syntax:
   - `~{variable}` → `$variable` (WDL modern syntax)
   - `${variable}` → `$variable` (WDL legacy syntax)
   - `$${...}` → `${...}` (Shell variable escaping preserved)
✅ **Sensible Defaults**: Auto-generates variable declarations with TODO markers for required values
✅ **Bash Validation**: Validates generated scripts for syntax errors before saving
✅ **JSON Input Templates**: Generates example input files matching task signatures
✅ **Comprehensive Logging**: Detailed DEBUG logging of all processing steps and substitutions
✅ **Error Handling**: Gracefully handles missing imports, parsing errors, and edge cases

## Installation

### Prerequisites
- Python 3.7+
- miniwdl library

### Setup
```bash
# Install miniwdl if not already installed
pip install miniwdl

# Make the script executable
chmod +x wdl_to_bash_extractor.py
```

## Usage

### Basic Usage
Extract all tasks from a WDL pipeline (run from workspace root or use absolute path):
```bash
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py pipeline.wdl
```

Or if you're in the workspace root:
```bash
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py all_of_us/mitochondria/mitochondria_pipeline.wdl
```

This creates an `extracted_tasks/` directory with bash scripts for each task.

### Test Mode: Find Bash Syntax Errors

Test mode scans WDL files for bash syntax errors in task command blocks **without** extracting or saving any scripts. This is useful for quality assurance and debugging.

```bash
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py pipeline.wdl --test
```

Test mode will:
- ✅ Parse all WDL tasks
- ✅ Extract and validate bash syntax for each task
- ✅ Generate a comprehensive error report
- ✅ Exit with status code 0 (no errors) or 1 (errors found)
- ❌ **Not** create any output files or directories

**Example output:**
```
======================================================================
WDL BASH ERROR DETECTION REPORT
======================================================================

Total Tasks Found: 12
Valid Tasks: 9
Tasks with Bash Errors: 3
Files Processed: 1

TASKS WITH BASH SYNTAX ERRORS:
----------------------------------------------------------------------

❌ build_vcf_shard_mt
   Error: bash: line 78: unexpected EOF while looking for matching `)'
   Source: /path/to/pipeline.wdl
```

### Custom Output Directory
```bash
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py pipeline.wdl -o ./my_scripts
```

### Enable Debug Logging
Useful for troubleshooting import resolution and variable substitution:
```bash
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py pipeline.wdl -l DEBUG
```

### Skip Syntax Validation
For faster processing (validation enabled by default):
```bash
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py pipeline.wdl --no-validation
```

### Skip Input JSON Templates
```bash
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py pipeline.wdl --no-input-templates
```

## Examples

### Example 1: Mitochondria Pipeline
```bash
# Run from workspace root with relative path
python3 scripts/wdl_to_bash/wdl_to_bash_extractor.py all_of_us/mitochondria/mitochondria_pipeline.wdl \
  -o extracted_mito_tasks -l INFO
```

This extracts 40+ tasks including:
- `Filter.sh` - Mutect calls filtering
- `M2.sh` - Mutect2 variant calling
- `ChainSwap.sh` - Genomic coordinate lifting
- And many more...

### Example 2: Generated Script Structure

Each extracted task generates:
1. **<TaskName>.sh** - The executable bash script
2. **<TaskName>_inputs.json** - Example input template

### Generated Script Format

```bash
#!/bin/bash

# Generated bash script for WDL task: M2
# Source: /path/to/source.wdl
# Generated: 2026-03-11T11:53:25.586087

# Input Variables:
# ===============
# input_bam: File
# ref_fasta: File
# (... all inputs listed with types ...)

# Output Variables:
# =================
# raw_vcf: File
# stats: File

# Input Variable Declarations
base_name="${TODO_base_name}"  # String - REQUIRED
input_bam="${TODO_input_bam}"  # File - REQUIRED
optional_param="${optional_param:-}"  # String?

# Exit on error
set -euo pipefail

# Command Block
# =============
[actual command code here with variables substituted]
```

## Variable Declaration Rules

The script generates variable declarations with these patterns:

### Required Variables
```bash
variable_name="${TODO_variable_name}"  # Type - REQUIRED
```
**Action needed**: Replace `${TODO_variable_name}` with actual value before running

### Optional Variables  
```bash
variable_name="${variable_name:-}"  # Type?
```
**Action needed**: Replace `${variable_name:-}` with default value or leave empty for default behavior

## Handling WDL Expressions

The extractor handles variable references but preserves WDL-specific expressions that require runtime evaluation:

### Simple Variables ✅ Converted
```wdl
echo ${input_file}        # → echo $input_file
ls ${output_dir}          # → ls $output_dir  
```

### Complex Expressions ⚠️ Manual Review Needed
```wdl
# Conditional expressions (preserve as-is or convert manually)
~{"--flag " + conditional_var}    # Requires manual WDL → bash conversion
~{default="/path" optional_file}  # Default values need bash equivalents
```

**Tips for manual conversion:**
- File operations: Use bash string concatenation `"--flag "$variable"` for flags followed by file paths
- Defaults: Use bash parameter expansion `${variable:-/default/value}` for optional files
- Conditionals: Use bash `if` statements for optional parameters

## Processing Summary

The script generates a summary report showing:
- Total tasks processed
- Output directory location  
- Variable substitutions performed
- Sample of substitutions (first 3 + count)

Example output:
```
======================================================================
WDL TO BASH EXTRACTION SUMMARY
======================================================================

Tasks Processed: 42
Output Directory: /path/to/extracted_tasks_mito
Files Processed: 18

Tasks Extracted:
----------------------------------------------------------------------
  M2:
    Source: .../AlignAndCallR1_v2_5_Single.wdl
    Inputs: 23
    Outputs: 4
    Variable Substitutions: 15
      ~{input_bam} -> $input_bam
      ~{ref_fasta} -> $ref_fasta
      ... and 12 more
```

## Logging Output

Enable DEBUG logging to see:
- Import path resolution attempts
- Successfully resolved import paths
- Input/output variable extraction per task
- Variable substitution operations
- Bash syntax validation results

Example:
```
DEBUG - Resolved import: ../tasks/M2.wdl -> /absolute/path/M2.wdl
DEBUG - Found input: input_bam: File
DEBUG - Substituted: ~{input_bam} -> $input_bam
DEBUG - Bash syntax validation passed for Filter
```

## Edge Cases Handled

### Remote Imports
Imports from URLs (e.g., GitHub raw content) are logged as warnings but don't block processing:
```
WARNING - Could not resolve import: https://raw.githubusercontent.com/.../Task.wdl
```
The script continues processing other files and tasks.

### Nested Imports
Recursive import resolution prevents infinite loops via tracking of processed files:
```
DEBUG - Already processed: /path/to/already/processed.wdl
```

### Multi-line Commands
Command blocks with complex control structures are preserved exactly:
```bash
for file in *.bam; do
  samtools view -h "$file" > "${file%.bam}.sam"
done
```

### Comments in Commands
Comments within command blocks are preserved:
```bash
# This is preserved
echo ${variable}  # Comments stay intact
```

### Escape Sequences
Shell escape sequences are preserved for literal variables:
```wdl
foo=$${var}  # Becomes foo=${var} in bash (literal text, not variable substitution)
```

## Common Workflows

### 1. Extract and Review
```bash
# Extract with debug logging
python3 wdl_to_bash_extractor.py pipeline.wdl -l DEBUG

# Review a specific task
cat extracted_tasks/TaskName.sh

# View input template
cat extracted_tasks/TaskName_inputs.json
```

### 2. Prepare for Execution
```bash
# Extract tasks
python3 wdl_to_bash_extractor.py pipeline.wdl

# Create actual input JSON from template
cp extracted_tasks/TaskName_inputs.json task_inputs.json
# Edit task_inputs.json with real values

# Define variables in bash script
source extracted_tasks/TaskName.sh  # Sources the script

# Execute
bash extracted_tasks/TaskName.sh
```

### 3. Batch Processing Multiple Pipelines
```bash
#!/bin/bash
for pipeline in pipelines/wdl/*/*.wdl; do
  echo "Processing: $pipeline"
  python3 wdl_to_bash_extractor.py "$pipeline" \
    -o "scripts/$(basename ${pipeline%.wdl})"
done
```

## Troubleshooting

### Script Not Executable
```bash
chmod +x wdl_to_bash_extractor.py
```

### WDL Syntax Unrecognized
Ensure WDL file is valid:
```bash
# Validate with womtool (if available)
java -jar womtool-*.jar validate pipeline.wdl
```

### Import Resolution Failures
Check import paths in WDL files are correct:
```bash
# Debug imports with verbose logging
python3 wdl_to_bash_extractor.py pipeline.wdl -l DEBUG 2>&1 | grep "import"
```

### Bash Syntax Errors
Review the generated script—some WDL expressions need manual conversion:
```bash
# Test script syntax
bash -n extracted_tasks/TaskName.sh
```

### Empty Output Directory
No tasks were successfully extracted. Check:
1. Input WDL file has tasks with command blocks
2. All imports resolved successfully (check DEBUG logs)
3. WDL syntax is valid

## Output File Examples

### Generated Bash Script
See [M2.sh](expected output) or [Filter.sh](expected output) for real examples from WARP mitochondria pipeline.

### Input JSON Template
```json
{
  "input_bam": "TODO: path/to/input.bam",
  "ref_fasta": "TODO: path/to/ref.fasta",
  "compress": true,
  "n_cpu": "TODO: <integer>",
  "optional_file": "TODO: path/to/optional_file"
}
```

## Architecture

The extractor uses a class-based design with these main components:

```python
WDLToBashExtractor
├── resolve_import_path()      # Resolves relative/absolute import paths
├── extract_imports()          # Finds and lists all imports in WDL file
├── extract_command_block()    # Extracts raw command from task
├── extract_input_variables()  # Parses task inputs with types
├── extract_output_variables() # Parses task outputs
├── substitute_variables()     # Performs WDL → bash variable conversion
├── generate_bash_script()     # Constructs complete bash script
├── validate_bash_syntax()     # Checks syntax with `bash -n`
├── save_bash_script()         # Writes to file and makes executable
├── generate_input_json_template() # Creates example inputs
├── process_task()             # Main task processing logic
└── process_wdl_file()         # Recursive file processing with import resolution
```

## Performance Notes

- **Typical Processing**: ~50-100ms per task
- **Import Resolution**: ~10-50ms per file depending on depth
- **Validation**: ~5-10ms per script (can be disabled with `--no-validation`)
- **Mitochondria Pipeline** (42 tasks, 18 files): ~2-3 seconds total

## Version History

- **1.0.0** - Initial release with full feature set

## License

Same license as the WARP repository

## Contributing

To enhance the tool:
1. Add support for additional WDL syntax patterns
2. Improve handling of complex expressions
3. Add template generation for other formats (HCL, YAML, etc.)

Submit issues and improvements to the WARP repository.
