version 1.0

# Pipeline for checking zero metrics in Optimus output files
workflow OptimusZeroMetricQC {
  input {
    File gene_metrics
    File cell_metrics  
    File library_metrics
    
    # Optional skip lists for customizing which metrics to ignore
    Array[String]? gene_skip_metrics
    Array[String]? cell_skip_metrics
    Array[String]? library_skip_metrics
  }

  call CheckForZeroColumns {
    input:
      gene_metrics = gene_metrics,
      cell_metrics = cell_metrics,
      library_metrics = library_metrics,
      gene_skip_metrics = gene_skip_metrics,
      cell_skip_metrics = cell_skip_metrics,
      library_skip_metrics = library_skip_metrics
  }

  output {
    String result = CheckForZeroColumns.result
    String gene_check_result = CheckForZeroColumns.gene_check_result
    String cell_check_result = CheckForZeroColumns.cell_check_result
    String library_check_result = CheckForZeroColumns.library_check_result
  }
}

# Task to check for zero values in Optimus metrics files
task CheckForZeroColumns {
  input {
    File gene_metrics
    File cell_metrics  
    File library_metrics
    
    # Optional skip lists
    Array[String]? gene_skip_metrics
    Array[String]? cell_skip_metrics
    Array[String]? library_skip_metrics
  }

  command <<<
    set -e
    pip3 install pandas
    python3 <<CODE
import pandas as pd
import sys

def check_zero_columns(file_path, file_name, skip_metrics=None):
    """Check for columns that are all zeros in a CSV/TSV file"""
    if skip_metrics is None:
        skip_metrics = [
            "reads_mapped_exonic_as"
        ]
    
    try:
        # Try reading as CSV first, then TSV
        try:
            df = pd.read_csv(file_path)
        except:
            df = pd.read_csv(file_path, sep='\t')
        
        # Find numeric columns that are all zeros
        zero_columns = []
        skipped_columns = []
        for col in df.columns:
            if pd.api.types.is_numeric_dtype(df[col]):
                if (df[col] == 0).all():
                    if col in skip_metrics:
                        skipped_columns.append(col)
                    else:
                        zero_columns.append(col)
        
        if skipped_columns:
            print(f"INFO: {file_name} skipped zero-value columns: {', '.join(skipped_columns)}")
        
        if zero_columns:
            print(f"ERROR: {file_name} has columns with all zero values: {', '.join(zero_columns)}")
            return False
        else:
            print(f"PASS: {file_name} has no problematic columns with all zero values")
            return True
            
    except Exception as e:
        print(f"ERROR: Could not process {file_name}: {str(e)}")
        return False

def check_zero_values_in_library_metrics(file_path, file_name, skip_metrics=None):
    """Check for metrics with zero values in a key-value format file"""
    if skip_metrics is None:
        skip_metrics = [
            "reads_mapped_antisense_to_gene",
            "reads_mapped_confidently_to_intronic_regions", 
            "percent_intronic_reads"
        ]
    
    try:
        zero_metrics = []
        skipped_metrics = []
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Parse key-value pairs (assuming format like "metric_name: value" or "metric_name\tvalue" or CSV)
                if ':' in line:
                    key, value = line.split(':', 1)
                elif '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        key, value = parts[0], parts[1]
                    else:
                        continue
                elif ',' in line:
                    # Handle CSV format (comma-separated)
                    parts = line.split(',')
                    if len(parts) >= 2:
                        key, value = parts[0], parts[1]
                    else:
                        continue
                else:
                    continue
                
                key = key.strip()
                value = value.strip()
                
                # Try to convert value to float and check if it's zero
                try:
                    numeric_value = float(value)
                    if numeric_value == 0.0:
                        if key in skip_metrics:
                            skipped_metrics.append(key)
                        else:
                            zero_metrics.append(key)
                except ValueError:
                    # Skip non-numeric values
                    continue
        
        if skipped_metrics:
            print(f"INFO: {file_name} skipped zero-value metrics: {', '.join(skipped_metrics)}")
        
        if zero_metrics:
            print(f"ERROR: {file_name} has metrics with zero values: {', '.join(zero_metrics)}")
            return False
        else:
            print(f"PASS: {file_name} has no problematic metrics with zero values")
            return True
            
    except Exception as e:
        print(f"ERROR: Could not process {file_name}: {str(e)}")
        return False

# Parse skip lists from WDL inputs
gene_skip = ~{if defined(gene_skip_metrics) then write_lines(select_first([gene_skip_metrics])) else ""}
cell_skip = ~{if defined(cell_skip_metrics) then write_lines(select_first([cell_skip_metrics])) else ""}
library_skip = ~{if defined(library_skip_metrics) then write_lines(select_first([library_skip_metrics])) else ""}

# Convert skip lists to Python lists
def parse_skip_list(skip_file_path):
    if not skip_file_path:
        return None
    try:
        with open(skip_file_path, 'r') as f:
            return [line.strip() for line in f if line.strip()]
    except:
        return None

gene_skip_list = parse_skip_list("${gene_skip}") if "${gene_skip}" else None
cell_skip_list = parse_skip_list("${cell_skip}") if "${cell_skip}" else None  
library_skip_list = parse_skip_list("${library_skip}") if "${library_skip}" else None

# Check each metrics file
print("=== OPTIMUS ZERO METRIC QC RESULTS ===")

gene_passed = check_zero_columns("~{gene_metrics}", "Gene Metrics", gene_skip_list)
cell_passed = check_zero_columns("~{cell_metrics}", "Cell Metrics", cell_skip_list) 
library_passed = check_zero_values_in_library_metrics("~{library_metrics}", "Library Metrics", library_skip_list)

all_passed = gene_passed and cell_passed and library_passed

# Write individual results to files
with open("gene_check_result.txt", "w") as f:
    f.write("PASS" if gene_passed else "FAIL")
    
with open("cell_check_result.txt", "w") as f:
    f.write("PASS" if cell_passed else "FAIL")
    
with open("library_check_result.txt", "w") as f:
    f.write("PASS" if library_passed else "FAIL")

print("=== OVERALL SUMMARY ===")
if not all_passed:
    print("OVERALL: One or more metrics files contain problematic zero values")
    sys.exit(1)
else:
    print("OVERALL: All metrics files passed zero value checks")

CODE
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/warp-tools/python-numpy-pandas:1.0.0-2.2.3-1.25.2"
    memory: "16 GiB"
    disks: "local-disk 10 HDD"
  }

  output {
    String result = "completed"
    String gene_check_result = read_string("gene_check_result.txt")
    String cell_check_result = read_string("cell_check_result.txt") 
    String library_check_result = read_string("library_check_result.txt")
  }
}
