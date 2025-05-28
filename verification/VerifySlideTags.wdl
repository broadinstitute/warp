version 1.0

import "../verification/VerifyTasks.wdl" as VerifyTasks

workflow VerifySlideTags {

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

        # spatial and positioning outputs 
        File test_spatial_output_h5
        File truth_spatial_output_h5

        File test_seurat_qs 
        File truth_seurat_qs 

        File test_coords_csv
        File truth_coords_csv

        File test_coords2_csv
        File truth_coords2_csv

        File test_intermediates_file
        File truth_intermediates_file

        Boolean? done
  }
    
    call VerifyTasks.CompareBams as CompareOptimusBams {
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
    
    call VerifyTasks.CompareH5Files as CompareSpatialOutputH5 {
        input:
            test_h5  = test_spatial_output_h5,
            truth_h5 = truth_spatial_output_h5
    }

    call compare_slidetags_csv as CompareCSV {
	      input:
              test_csv  = test_coords_csv,
              truth_csv = truth_coords_csv
    }
    
    call compare_slidetags_csv2 as CompareCSV2 {
	      input:
              test_csv  = test_coords2_csv,
              truth_csv = truth_coords2_csv
    }

    call CompareSlideTagsTarballContents as CompareTarballContents {
	      input:
		        test_tar  = test_intermediates_file,
		        truth_tar = truth_intermediates_file
    }

}

task CompareSlideTagsTarballContents {
  input {
    File test_tar
    File truth_tar
  }

  Float file_size = size(test_tar, "GiB") + size(truth_tar, "GiB")
  Int disk_size = ceil(file_size * 4) + 20

  command <<<
    set -e

    mkdir test_dir truth_dir

    # Extract tarballs
    tar -xzf ~{test_tar} -C test_dir
    tar -xzf ~{truth_tar} -C truth_dir

    # Prepare log file
    touch comparison_errors.log

    # Compare matrix.csv.gz
    echo "Comparing matrix.csv.gz..."
    gunzip -c test_dir/matrix.csv.gz | sort > sorted_test_matrix.txt
    gunzip -c truth_dir/matrix.csv.gz | sort > sorted_truth_matrix.txt
    if ! diff sorted_test_matrix.txt sorted_truth_matrix.txt > /dev/null; then
        echo "matrix.csv.gz files differ." >> comparison_errors.log
    else
        echo "matrix.csv.gz files are identical."
    fi

    # Compare cb_whitelist.txt
    echo "Comparing cb_whitelist.txt..."
    sort test_dir/cb_whitelist.txt > sorted_test_cb.txt
    sort truth_dir/cb_whitelist.txt > sorted_truth_cb.txt
    if ! diff sorted_test_cb.txt sorted_truth_cb.txt > /dev/null; then
        echo "cb_whitelist.txt files differ." >> comparison_errors.log
    else
        echo "cb_whitelist.txt files are identical."
    fi

    # Final check for errors
    if [ -s comparison_errors.log ]; then
        sed -i '1iNOTE: These differences may be due to non-deterministic behavior in the pipeline, but should still be investigated.' comparison_errors.log
        echo "Comparison failed. See comparison_errors.log for details:"
        cat comparison_errors.log
        exit 1
    else
        echo "All comparisons succeeded."
    fi
  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk ~{disk_size} HDD"
    memory: "20 GiB"
    preemptible: 3
  }
}


task compare_slidetags_csv {
  input {
    File truth_csv
    File test_csv
  }

  command <<<
python3 <<CODE
import csv
import sys

truth_file = "~{truth_csv}"
test_file = "~{test_csv}"

columns_to_compare = ["cb", "umi", "beads", "max"]
mismatches = []

with open(truth_file, newline='') as f1, open(test_file, newline='') as f2:
    reader1 = csv.reader(f1)
    reader2 = csv.reader(f2)

    header1 = next(reader1)
    header2 = next(reader2)

    if header1 != header2:
        print("ERROR: CSV headers do not match.")
        print("Truth header:", header1)
        print("Test header:", header2)
        sys.exit(1)

    col_indices = {name: header1.index(name) for name in columns_to_compare}

    row_num = 1
    for row1, row2 in zip(reader1, reader2):
        row_num += 1
        for colname in columns_to_compare:
            idx = col_indices[colname]
            val1 = row1[idx]
            val2 = row2[idx]

            print(f"Checking {colname}:")
            print(f"  Truth Row {row_num}: {val1}")
            print(f"  Test Row {row_num}: {val2}")

            if val1 != val2:
                print(f"  --> Mismatch: {val1} != {val2}")
                mismatches.append((row_num, colname, val1, val2))

if mismatches:
    print("\nSummary of mismatches (exact match expected):")
    for row_num, col, v1, v2 in mismatches:
        print(f"- Row {row_num}, Column {col}: {v1} != {v2}")
    print("Comparison failed.")
    sys.exit(1)
else:
    print("Files match exactly for cb, umi, beads, and max columns.")
CODE
  >>>

  runtime {
    docker: "python:3.9"
  }
}

task compare_slidetags_csv2 {
  input {
    File truth_csv
    File test_csv
  }

  command <<<
python3 <<CODE
import csv
import sys

truth_file = "~{truth_csv}"
test_file = "~{test_csv}"

columns_to_compare = ["cb", "umi", "beads"]
mismatches = []

with open(truth_file, newline='') as f1, open(test_file, newline='') as f2:
    reader1 = csv.reader(f1)
    reader2 = csv.reader(f2)

    header1 = next(reader1)
    header2 = next(reader2)

    if header1 != header2:
        print("ERROR: CSV headers do not match.")
        print("Truth header:", header1)
        print("Test header:", header2)
        sys.exit(1)

    col_indices = {name: header1.index(name) for name in columns_to_compare}

    row_num = 1
    for row1, row2 in zip(reader1, reader2):
        row_num += 1
        for colname in columns_to_compare:
            idx = col_indices[colname]
            val1 = row1[idx]
            val2 = row2[idx]

            print(f"Checking {colname}:")
            print(f"  Truth Row {row_num}: {val1}")
            print(f"  Test Row {row_num}: {val2}")

            if val1 != val2:
                print(f"  --> Mismatch: {val1} != {val2}")
                mismatches.append((row_num, colname, val1, val2))

if mismatches:
    print("\nSummary of mismatches (exact match expected):")
    for row_num, col, v1, v2 in mismatches:
        print(f"- Row {row_num}, Column {col}: {v1} != {v2}")
    print("Comparison failed.")
    sys.exit(1)
else:
    print("Files match exactly for cb, umi, and beads columns.")
CODE
  >>>

  runtime {
    docker: "python:3.9"
  }
}


