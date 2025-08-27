version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

# task to warn if optimus metrics are all zero
task CheckForZeroColumns {
  input {
    File gene_metrics
    File cell_metrics  
    File library_metrics
  }

  command <<<
    set -e
    pip3 install pandas
    python3 <<CODE
import pandas as pd
import sys

def check_zero_columns(file_path, file_name):
    """Check for columns that are all zeros in a CSV/TSV file"""
    try:
        # Try reading as CSV first, then TSV
        try:
            df = pd.read_csv(file_path)
        except:
            df = pd.read_csv(file_path, sep='\t')
        
        # Find numeric columns that are all zeros
        zero_columns = []
        for col in df.columns:
            if pd.api.types.is_numeric_dtype(df[col]):
                if (df[col] == 0).all():
                    zero_columns.append(col)
        
        if zero_columns:
            print(f"ERROR: {file_name} has columns with all zero values: {', '.join(zero_columns)}")
            return False
        else:
            print(f"PASS: {file_name} has no columns with all zero values")
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

# Check each metrics file
all_passed = True

all_passed &= check_zero_columns("~{gene_metrics}", "Gene Metrics")
all_passed &= check_zero_columns("~{cell_metrics}", "Cell Metrics") 
all_passed &= check_zero_values_in_library_metrics("~{library_metrics}", "Library Metrics")

if not all_passed:
    print("OVERALL: One or more metrics files contain columns with all zero values")
    sys.exit(1)
else:
    print("OVERALL: All metrics files passed zero column checks")

CODE
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/warp-tools/python-numpy-pandas:1.0.0-2.2.3-1.25.2"
    memory: "16 GiB"
    disks: "local-disk 10 HDD"
  }

  output {
    String result = "completed"
  }
}

workflow VerifyOptimus {
  
  input {
    File test_h5ad
    File truth_h5ad

    File test_bam
    File truth_bam

    File test_gene_metrics
    File truth_gene_metrics
    
    File test_cell_metrics
    File truth_cell_metrics

    File test_library_metrics
    File truth_library_metrics

  Boolean check_zero_metrics = false
    Boolean? done
  }

  call VerifyTasks.CompareBams as CompareBams {
    input:
      test_bam       = test_bam,
      truth_bam      = truth_bam,
      lenient_header = true
  }

  call VerifyTasks.CompareCompressedTextFiles as CompareGeneMetrics {
    input:
      test_zip  = test_gene_metrics,
      truth_zip = truth_gene_metrics
  }

  call VerifyTasks.CompareCompressedTextFiles as CompareCellMetrics {
    input:
      test_zip  = test_cell_metrics,
      truth_zip = truth_cell_metrics
  }
  
  call VerifyTasks.CompareH5adFilesGEX as CompareH5adFilesOptimus {
        input:
            test_h5ad  = test_h5ad,
            truth_h5ad = truth_h5ad
    }

  call VerifyTasks.CompareLibraryFiles as CompareLibraryMetrics {
          input:
            test_text_file = test_library_metrics,
            truth_text_file = truth_library_metrics
  }

  if (check_zero_metrics) {
    call CheckForZeroColumns as CheckZeroColumns {
      input:
        gene_metrics = test_gene_metrics,
        cell_metrics = test_cell_metrics,
        library_metrics = test_library_metrics
    }
  }

  output {
    String? zero_check_result = CheckZeroColumns.result
  }
}