# WDL to Bash Extractor - Advanced Features & Edge Cases

## Advanced Features

### 1. Recursive Import Resolution

The extractor automatically follows `import` statements across multiple files with cycle detection.

**How it works:**
- Finds all `import "path/to/file.wdl"` statements
- Resolves relative paths from the importing file's directory
- Processes imports depth-first (imported files first)
- Prevents infinite loops by tracking processed files

**Example:**
```
mitochondria_pipeline.wdl
├── imports AlignAndCallR1_v2_5_Single.wdl
│   ├── imports MongoTasks_v2_5_Single.wdl
│   │   └── (process tasks here)
│   └── (process tasks)
├── imports ProduceSelfReferenceFiles_v2_5_Single.wdl
│   └── (process tasks)
└── (process tasks)
```

**Result:** Tasks from all files available in single output directory

### 2. Type-Aware Variable Generation

The extractor understands WDL types and generates appropriate bash variable patterns.

**Type Mapping:**
```
WDL Type              | Bash Pattern                 | Example
String                | ${TODO_varname}              | "${TODO_name}"
File                  | ${TODO_varname}              | "${TODO_input_file}"
Int                   | ${TODO_varname}              | "${TODO_threads}"
Float                 | ${TODO_varname}              | "${TODO_threshold}"
Boolean               | ${TODO_varname}              | "${TODO_run_qc}"
String?               | ${varname:-}                 | "${optional_str:-}"
Array[String]         | ${varname:-}                 | "${array_var:-}"
```

**Optional Types:** Automatically detected and marked with `${var:-}` default pattern

### 3. Comprehensive Variable Substitution

Handles both modern (`~{...}`) and legacy (`${...}`) WDL variable syntax.

**Substitution Process:**
1. Protects shell escapes (`$${...}`) by temporary replacement
2. Replaces WDL syntax for all known input variables
3. Restores shell escapes as bash syntax (`${...}`)

**Example:**
```wdl
# Input: WDL command block
echo ~{input_name}
cat ${data_file}
export HOME=$${HOME}  # Shell variable

# Output: Bash command block
echo $input_name
cat $data_file
export HOME=${HOME}  # Shell variable (preserved)
```

### 4. Bash Syntax Validation

Each generated script is validated using `bash -n` before saving.

**What gets validated:**
- Correct bracket matching
- Valid bash syntax
- Function declarations
- Control structure syntax (if/for/while)

**Validation errors are:**
- Logged with details
- Prevent script from being saved
- Can be skipped with `--no-validation` flag

### 5. JSON Input Templates

Auto-generates input files with appropriate defaults based on types.

**Default Patterns:**
```json
{
  "file_variable": "TODO: path/to/file",
  "string_variable": "TODO: string_value",
  "integer_variable": "TODO: <integer>",
  "float_variable": "TODO: <float>",
  "bool_variable": true,
  "array_variable": ["TODO: array_element"]
}
```

## Edge Cases Handled

### Case 1: Remote Imports (URLs)

**Scenario:** WDL imports from GitHub or other remote sources
```wdl
import "https://raw.githubusercontent.com/org/repo/main/Task.wdl" as Task
```

**Handling:**
- Logged as WARNING (import could not be resolved)
- Processing continues with other files
- Tasks from that import skipped
- No blocker to extraction success

**Output:**
```
WARNING - Could not resolve import: https://raw.githubusercontent.com/...
```

### Case 2: Missing Files

**Scenario:** Imported file doesn't exist locally
```wdl
import "../missing/Task.wdl" as Task
```

**Handling:**
- Path resolution fails
- Logged as WARNING
- Extraction continues with available files
- User notified in summary report

**Debug:**
```bash
python3 wdl_to_bash_extractor.py pipeline.wdl -l DEBUG 2>&1 | grep -i "import\|resolv"
```

### Case 3: Circular Imports

**Scenario:** Files import each other (A→B→A)
```wdl
# A.wdl
import "B.wdl" as B

# B.wdl
import "A.wdl" as A
```

**Handling:**
- Cycle detection via file path tracking
- Second occurrence logged as DEBUG skip
- No infinite loop or duplicate processing

**Output:**
```
DEBUG - Already processed: /path/to/A.wdl
```

### Case 4: Tasks Without Command Blocks

**Scenario:** WDL has declarations but no executable command
```wdl
task DeclarativeTask {
  input {
    String value
  }
  # No command block - this is a compute task or pure declaration
}
```

**Handling:**
- Skipped with WARNING message
- Not counted in final extraction count
- Other tasks continue normally

**Output:**
```
WARNING - Task DeclarativeTask has no command block, skipping
```

## Complex WDL Features

### Feature 1: WDL Expressions in Commands

Some WDL features require manual conversion (extractor logs warning).

**Type: Conditional Expressions**
```wdl
# Original WDL
command <<<
  gatk FilterMutectCalls \
    ~{"--min-allele-fraction " + vaf_threshold} \
    ~{"--contamination " + contamination}
>>>

# Extracted bash (needs manual review + conversion)
gatk FilterMutectCalls \
  ~{"--min-allele-fraction " + vaf_threshold} \
  ~{"--contamination " + contamination}

# Manual bash conversion
gatk FilterMutectCalls \
  ${vaf_threshold:+--min-allele-fraction $vaf_threshold} \
  ${contamination:+--contamination $contamination}
```

**Type: Default Values**
```wdl
# Original WDL
command <<<
  export PATH=~{default="/usr/bin" custom_path}:$PATH
>>>

# Manual bash conversion
export PATH="${custom_path:-/usr/bin}:$PATH"
```

### Feature 2: Runtime Sections

WDL runtime sections are extracted as comments (not executable in bash).

```wdl
# Original WDL
task AlignTask {
  command { ... }
  runtime {
    docker: "broadinstitute/gatk:4.2.0"
    cpu: 8
    memory: "16 GB"
    disk: "100 GB"
  }
}

# Extracted bash (runtime info in comments/header only)
# Runtime Requirements:
# - Docker: broadinstitute/gatk:4.2.0
# - CPU: 8 cores
# - Memory: 16 GB
# - Disk: 100 GB
```

**To enforce runtime requirements in bash, create Docker wrapper:**
```bash
#!/bin/bash
docker run --rm \
  -v "$PWD:$PWD" \
  -w "$PWD" \
  --cpus=8 \
  --memory=16g \
  broadinstitute/gatk:4.2.0 \
  bash extracted_tasks/TaskName.sh
```

### Feature 3: WDL Variable Interpolation

Some interpolation patterns need special handling.

**Simple Interpolation (Handled):**
```wdl
echo "Files: ${file1}, ${file2}"  # → echo "Files: $file1, $file2"
```

**String Concatenation (Manual):**
```wdl
# Original WDL
echo "Output: ${output_dir}/${sample_name}.vcf"

# Extracted (needs bash fixing):
echo "Output: $output_dir$sample_name.vcf"

# Should be (add slashes):
echo "Output: $output_dir/$sample_name.vcf"
```

## Performance Optimization

### For Very Large Pipelines

**Strategy 1: Skip Expensive Operations**
```bash
# Disable bash validation (5-10ms per script)
python3 wdl_to_bash_extractor.py massive_pipeline.wdl --no-validation

# Disable JSON generation (2-3ms per task)
python3 wdl_to_bash_extractor.py massive_pipeline.wdl --no-input-templates

# Combine both
python3 wdl_to_bash_extractor.py massive_pipeline.wdl \
  --no-validation --no-input-templates
```

**Sample timings (mitochondria pipeline with 42 tasks):**
- Full extraction with validation: ~2-3 seconds
- Without validation: ~1-2 seconds
- Without JSON: ~1-2 seconds
- Both disabled: ~0.5-1 second

**Strategy 2: Parallel Processing**
```bash
#!/bin/bash
# Extract multiple pipelines in parallel

pipelines=(
  "pipelines/wdl/dna_seq/GermlineSingleSample.wdl"
  "pipelines/wdl/rna_seq/Optimus.wdl"
  "all_of_us/mitochondria/mitochondria_pipeline.wdl"
)

export -f extract_pipeline
for pipeline in "${pipelines[@]}"; do
  python3 wdl_to_bash_extractor.py "$pipeline" \
    -o "extracted/$(basename ${pipeline%.wdl})" &
done

wait  # Wait for all background jobs
```

## Troubleshooting Advanced Issues

### Issue: Variable Substitution Incomplete

**Symptom:** Some `~{...}` or `${...}` remain in script

**Cause:** Variables not in input section (computed/derived variables)

**Solution:**
```bash
# Check for non-input variables in command
grep -E '\$\{[A-Za-z_]' extracted_tasks/TaskName.sh | grep -v "^export"

# Manually add computed variables
# Example: if command uses ${output_filename} but not in inputs
output_filename="${sample_name}.vcf"
```

### Issue: Import Resolution Failing

**Symptom:** 
```
WARNING - Could not resolve import: ../relative/path/Task.wdl
```

**Diagnosis:**
```bash
# Verify relative path from importing file
cd path/to/importing/directory
ls -la ../relative/path/Task.wdl

# Check for symbolic links
readlink -f ../relative/path/Task.wdl
```

**Fix:**
- Ensure WDL uses correct relative paths
- Make sure working directory when running extractor is correct
- Try absolute paths if relative resolution fails

### Issue: Circular Dependency in Imports

**Symptom:** Error during parsing or infinite processing

**Solution:**
```bash
# Extract with debug to see circular import
python3 wdl_to_bash_extractor.py pipeline.wdl -l DEBUG 2>&1 | \
  grep -E "Already processed|^Processing"

# Manually break cycle by restructuring imports
# Move shared tasks to separate file that both import from
```

### Issue: Bash Validation Failures

**Symptom:**
```
ERROR - Bash syntax error in TaskName: ...
Skipping TaskName due to bash syntax error
```

**Cause:** WDL-specific syntax remaining in script

**Common patterns:**
```bash
# WDL conditional - needs manual conversion
~{if condition then "true" else "false"}

# WDL ternary
condition ? "yes" : "no"

# WDL string interpolation with operations
"${file}" + ".backup"
```

**Fix:**
```bash
# Manually convert after extraction
nano extracted_tasks/TaskName.sh

# Convert WDL conditionals to bash
# ~{condition} → ${condition:+value}

# Test after manual edits
bash -n extracted_tasks/TaskName.sh
```

## Advanced Logging

### Full Diagnostic Logging

```bash
# Capture all log output
python3 wdl_to_bash_extractor.py pipeline.wdl -l DEBUG 2>&1 | tee extraction.log

# Analyze logs
echo "=== Import Resolution ==="
grep -i "import\|resolv" extraction.log

echo "=== Variable Processing ==="
grep -i "input\|output\|substitut" extraction.log

echo "=== Errors and Warnings ==="
grep -iE "error|warning" extraction.log
```

### Custom Logging Configuration

To modify logging behavior, edit the script's logging configuration:

```python
# In wdl_to_bash_extractor.py - line ~51
logging.basicConfig(
    level=getattr(logging, log_level),
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    # Add file handler:
    handlers=[
        logging.FileHandler('extraction_debug.log'),
        logging.StreamHandler()
    ]
)
```

## Integration Examples

### Integration with Nextflow

```groovy
// Nextflow process using extracted bash script
process runExtractedTask {
    input:
        file(bam_file)
        val(sample_id)
    output:
        file("*.vcf")
    script:
        """
        export input_bam=${bam_file}
        export sample_name=${sample_id}
        bash ${projectDir}/extracted_tasks/M2.sh
        """
}
```

### Integration with Snakemake

```python
# Snakemake rule using extracted bash script
rule run_extracted_task:
    input:
        bam="samples/{sample}.bam"
    output:
        vcf="variants/{sample}.vcf"
    shell:
        """
        export input_bam={input.bam}
        export output_file={output.vcf}
        bash extracted_tasks/TaskName.sh
        """
```

### Integration with Make

```makefile
# GNU Make using extracted scripts
extracted_tasks: 
	python3 wdl_to_bash_extractor.py pipeline.wdl -o $@

run_m2: extracted_tasks
	export input_bam=input.bam; \
	export ref_fasta=ref.fasta; \
	bash extracted_tasks/M2.sh

run_filter: extracted_tasks
	bash extracted_tasks/Filter.sh
```

## Best Practices

1. **Always validate after extraction**
   ```bash
   for script in extracted_tasks/*.sh; do
     bash -n "$script" || echo "ERROR in $script"
   done
   ```

2. **Use absolute paths for inputs**
   ```bash
   # Good
   input_bam="/data/samples/sample.bam"
   
   # Problematic (relative paths)
   input_bam="sample.bam"
   ```

3. **Document variable sources**
   ```bash
   # In your wrapper script, document where each variable comes from
   input_bam="${DATA_DIR}/samples/${SAMPLE_ID}.bam"  # From argument
   ref_fasta="/ref/genomes/hg38/fasta"               # Shared resource
   ```

4. **Test with sample data first**
   ```bash
   # Create minimal test inputs
   bash extracted_tasks/TaskName.sh 2>&1 | head -20
   ```

5. **Version control extracted scripts**
   ```bash
   git add extracted_tasks/
   git commit -m "Extract WDL tasks from pipeline v1.2.3"
   ```
