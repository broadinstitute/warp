version 1.0

task CompareVcfs {
  input {
    File file1
    File file2
    String patternForLinesToExcludeFromComparison = ""
  }

  command {
    set -eo pipefail

    if [ -z '~{patternForLinesToExcludeFromComparison}' ]; then
      diff <(gunzip -c -f ~{file1}) <(gunzip -c -f ~{file2})
    else
      echo "It's defined!"
      diff <(gunzip -c -f ~{file1} | grep -v '~{patternForLinesToExcludeFromComparison}') <(gunzip -c -f ~{file2} | grep -v '~{patternForLinesToExcludeFromComparison}')
    fi
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 70 HDD"
    memory: "32 GiB"
    preemptible: 3
  }
}

task CompareVcfsAllowingQualityDifferences {
  input {
    File file1
    File file2
  }

  command {
    exit_code=0

    cmp <(gunzip -c -f ~{file1} | grep -v '^##') <(gunzip -c -f ~{file2} | grep -v '^##')
    if [ $? -ne 0 ]; then
      exit_code=1
      echo "Error: VCF ~{file1} differs in content from ~{file2}" >&2
      cmp <(gunzip -c -f ~{file1} | grep -v '^##' | cut -f 1-5,7-) <(gunzip -c -f ~{file2} | grep -v '^##' | cut -f 1-5,7-)
      if [ $? -eq 0 ]; then
        echo "However they ONLY differ in the quality column" >&2
      fi
    fi

    exit $exit_code
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 50 HDD"
    memory: "3 GiB"
    preemptible: 3
  }
}

task CompareVCFsVerbosely {
  input {
    File actual
    File actual_index
    File expected
    File expected_index
    File ref_fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File ref_fasta_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    File ref_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    String extra_args = " --ignore-attribute VQSLOD --ignore-attribute AS_VQSLOD --ignore-filters "
    + "--ignore-attribute culprit --ignore-attribute AS_culprit --ignore-attribute AS_FilterStatus "
    + "--ignore-attribute ExcessHet --ignore-star-attributes --allow-nan-mismatch --ignore-attribute END"
    Boolean warn_on_error = true
  }

  command {
    gatk --java-options "-Xms2500m -Xmx2500m" VCFComparator -R ~{ref_fasta}  -V:actual ~{actual} -V:expected ~{expected} ~{extra_args} ~{if(warn_on_error) then "--warn-on-errors" else ""} --finish-before-failing
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    disks: "local-disk 50 HDD"
    memory: "3 GiB"
    preemptible: 3
  }
}

task CompareGtcs {
  input {
    File file1
    File file2
    File bead_pool_manifest_file
  }

  command {
    java -Xms4500m -Xmx4500m -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
      CompareGtcFiles \
      --INPUT ~{file1} \
      --INPUT ~{file2} \
      --BPM_FILE ~{bead_pool_manifest_file}
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk 10 HDD"
    memory: "5000 MiB"
    preemptible: 3
  }
}

task CompareTabix {
  input {
    File test_fragment_file
    File truth_fragment_file
  }
  command <<<
  exit_code=0

  a=$(gunzip -c ~{test_fragment_file} | md5sum | awk '{ print $1 }')
  b=$(gunzip -c ~{truth_fragment_file} | md5sum | awk '{ print $1 }')

  if [[ "$a" = "$b" ]]; then
    echo "The fragment files are equal"
  else
    echo "The fragment files md5sums do not match. Performing a line count:"

    test_lines=$(gunzip -c ~{test_fragment_file} | wc -l)
    truth_lines=$(gunzip -c ~{truth_fragment_file} | wc -l)

    echo "Test file has $test_lines lines"
    echo "Truth file has $truth_lines lines"

    diff_lines=$((test_lines - truth_lines))
    abs_diff_lines=${diff_lines#-}

    if [[ $abs_diff_lines -gt 100 ]]; then
      echo "Line count difference greater than 100 lines. The line count difference is $abs_diff_lines lines. Task failed."
      exit_code=1
    fi
  fi

  exit $exit_code

  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/snapatac2:2.0.0"
    disks: "local-disk 100 HDD"
    memory: "50 GiB"
    preemptible: 3
  }   
}

task CompareTextFiles {
  input {
    Array[File] test_text_files
    Array[File] truth_text_files
  }

  command {
    exit_code=0

    test_files_length=~{length(test_text_files)}
    truth_files_length=~{length(truth_text_files)}
    if [ $test_files_length -ne $truth_files_length ]; then
      exit_code=1
      echo "Error: Different number of input files ($test_files_length vs. $truth_files_length).  This is really not OK"
    fi

    while read -r a && read -r b <&3;
    do
      echo "Sorting File $a and $b"
      sort $a > $a.sorted
      sort $b > $b.sorted

      echo "Calculating md5sums for $a and $b"
      md5_a=$(md5sum $a.sorted | cut -d ' ' -f1)
      md5_b=$(md5sum $b.sorted | cut -d ' ' -f1)

      if [ $md5_a = $md5_b ]; then
        echo "Files $a.sorted and $b.sorted have matching md5sums and are the same."
      else
        echo "Files $a.sorted and $b.sorted have different md5sums."
        diff $a.sorted $b.sorted >&2
        exit_code=1
      fi

    done < ~{write_lines(test_text_files)} 3<~{write_lines(truth_text_files)}

    echo "Exiting with code $exit_code"
    exit $exit_code
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 100 HDD"
    memory: "50 GiB"
    preemptible: 3
  }
}


task CompareAtacLibraryMetrics {
  input {
    Array[File] test_text_files
    Array[File] truth_text_files
  }

  command <<<
python3 <<CODE
import csv
import hashlib

# Define acceptable percentage-based thresholds for nondeterministic metrics
# Arrived at these thresholds by examining the differences between the test and truth files in our scientific tests
thresholds = {
    "sequenced_reads": 0.0000000066,
    "fraction_Q30_bases_in_read_1": 0.0000000054,
    "fraction of high-quality fragments in cells": 0.000000054,
    "fraction_of_transposition_events_in_peaks_in_cells": 0.00000037,
    "fraction_duplicates": 0.00000017,
    "fraction_confidently_mapped": 0.000000123,
    "fraction_unmapped": 0.0000016,
    "fraction_nonnuclear": 0.00000079,
    "fraction_fragment_in_nucleosome_free_region": 0.00000059,
    "fraction_fragment_flanking_single_nucleosome": 0.00000057,
    "tss_enrichment_score": 0.0000024,
    "fraction_of_high-quality_fragments_overlapping_tss": 0.00000025,
    "number_of_peaks": 0.0000074,
    "fraction_of_genome_in_peaks": 0.0000024,
    "fraction_of_high-quality_fragments_overlapping_peaks": 0.00000030
}

thresholds = {k.lower(): v for k, v in thresholds.items()}


def calculate_md5(file_path):
    """Calculates the MD5 checksum for a file."""
    print(f"Processing file: {file_path}")
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def compare_files(test_file, truth_file):
    """Compare two files by their md5sums and then line by line for metric comparison."""
    test_md5 = calculate_md5(test_file)
    truth_md5 = calculate_md5(truth_file)

    if test_md5 == truth_md5:
        print(f"Files {test_file} and {truth_file} have matching md5sums. No further checks needed.")
        return True
    else:
        print(f"Files {test_file} and {truth_file} have different md5sums. Proceeding with metric comparison.")
        return compare_metrics(test_file, truth_file)

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def compare_metrics(test_file, truth_file):
    exit_code = 0
    with open(test_file, newline='') as test_f, open(truth_file, newline='') as truth_f:
        test_reader = csv.reader(test_f)
        truth_reader = csv.reader(truth_f)
        for (test_row, truth_row) in zip(test_reader, truth_reader):
            metric_a, value_a = test_row
            metric_b, value_b = truth_row

            # Skip non-numeric values
            if not is_float(value_a) or not is_float(value_b):
                print(f"Skipping non-numeric metric: {metric_a} or {metric_b}")
                continue

            value_a, value_b = float(value_a), float(value_b)
            if metric_a != metric_b:
                print(f"Error: Metric names don't match for {metric_a} and {metric_b}")
                exit_code = 1
                continue
            # Check if the metric has a set threshold, otherwise default to 0.00
            threshold = thresholds.get(metric_a.lower(), 0.00)

            diff = abs(value_a - value_b)

            # Calculate the allowable difference based on the threshold
            allowable_diff = value_b * threshold
            if diff > allowable_diff:
                print(f"Error: Metric {metric_a} exceeds threshold. Test value: {value_a}, Truth value: {value_b}, Threshold: {threshold*100}%. The allowable difference is {allowable_diff} and the actual difference is {diff}.")
                exit_code = 1
            else:
                print(f"Metric {metric_a} is within the threshold.")
    return exit_code == 0


# Read and compare all files
test_files = ["~{sep=',' test_text_files}"]
truth_files = ["~{sep=',' truth_text_files}"]

if len(test_files) != len(truth_files):
    print(f"Error: Different number of input files ({len(test_files)} vs. {len(truth_files)}). This is really not OK")
    exit(1)

for test_file, truth_file in zip(test_files, truth_files):
    if not compare_files(test_file, truth_file):
        exit(1)

print("All files passed the comparison.")
CODE

  >>>

  runtime {
    docker: "python:3.9-slim"
    disks: "local-disk 100 HDD"
    memory: "50 GiB"
    preemptible: 3
  }
}



task CompareCrams {

  input {
    File test_cram
    File test_crai
    File truth_cram
    File truth_crai

    Int disk_size_gb = ceil((size(test_cram, "GiB") + size(truth_cram, "GiB"))) + 50
  }

  command {
    # get the offset of the first alignment
    test_offset="$(zcat ~{test_crai} | cut -f4 | head -n 1)"
    truth_offset="$(zcat ~{truth_crai} | cut -f4 | head -n 1)"

    # compare files with byte offset
    cmp -i "$test_offset:$truth_offset" ~{test_cram} ~{truth_cram}
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk " + disk_size_gb + " HDD"
    memory: "2 GiB"
    preemptible: 3
  }
}

task CompareCrais {

  input {
    File test_crai
    File truth_crai
  }

  command {
    # compare columns 1,2,3,5, and 6. Cannot compare column 4
    # because it is a byte offset number that may differ.
    cmp <(zcat ~{test_crai} | cut -f1,2,3,5,6) <(zcat ~{truth_crai} | cut -f1,2,3,5,6)
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
}

task CompareBams {

  input {
    File test_bam
    File truth_bam
    Boolean lenient_header = false
    Boolean lenient_low_mq = false
  }

  Float bam_size = size(test_bam, "GiB") + size(truth_bam, "GiB")
  Int disk_size = ceil(bam_size * 4) + 200
  Int memory_mb = 600000
  Int java_memory_size = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    set -e
    set -o pipefail

    truth_bam=~{truth_bam}
    test_bam=~{test_bam}

    # Get the sizes of the BAM files in bytes
    truth_size=$(stat -c %s ~{truth_bam})
    test_size=$(stat -c %s ~{test_bam})

    # Calculate the difference in bytes
    size_difference=$((truth_size - test_size))

    # Calculate the absolute value of the difference
    abs_size_difference=$((size_difference < 0 ? -size_difference : size_difference))

    # Compare the sizes and fail fast if the difference is greater than 200 * 1024 * 1024 bytes (200 MB)
    if [ "$abs_size_difference" -gt $((200 * 1024 * 1024)) ]; then
        echo "Skipping CompareSAMs as BAM file sizes differ by more than 200 MB. $truth_bam is $truth_size bytes and $test_bam is $test_size bytes. Exiting."
        exit 1
    else
        echo "WARNING: BAM file sizes differ by less than 200 MB. $truth_bam is $truth_size bytes and $test_bam is $test_size bytes. Proceeding to CompareSAMs:"

        java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
        CompareSAMs \
            ~{test_bam} \
            ~{truth_bam} \
            O=comparison.tsv \
            LENIENT_HEADER=~{lenient_header} \
            LENIENT_LOW_MQ_ALIGNMENT=~{lenient_low_mq} \
            MAX_RECORDS_IN_RAM=300000
    fi

  >>>


  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    cpu: 2
    memory: "${memory_mb} MiB"
    preemptible: 3
  }
}

task CompareCompressedTextFiles {

  input {
    File test_zip
    File truth_zip
  }

  Float file_size = size(test_zip, "GiB") + size(truth_zip, "GiB")
  Int disk_size = ceil(file_size * 4) + 20

  command {
    diff <(gunzip -c ~{test_zip} | sort) <(gunzip -c ~{truth_zip} | sort)

    if [ $? -eq 0 ]; then
        echo "Comparison succeeded: The files are identical."
    else
        echo "Comparison failed: The files differ."
        exit 1
    fi
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk " + disk_size + " HDD"
    memory: "20 GiB"
    preemptible: 3
  }

}

task CompareLooms {

  input {
    File truth_loom
    File test_loom
    Float delta_threshold = 0.05

    Int cpu = 3
    String docker = "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    Int disk_size_gb = ceil((size(truth_loom, "GiB") + size(test_loom, "GiB")) * 2) + 20
    Int memory_mb = ceil(size(truth_loom, "MiB") + size(test_loom, "MiB") * 4) + 20000
  }

  command <<<
  set -e
  pip3 install 'matplotlib<3.7' scanpy loompy numpy pandas

  python3 <<CODE
  import sys
  import scanpy
  import numpy as np

  test_loom_file="~{test_loom}"
  truth_loom_file = "~{truth_loom}"
  threshold = ~{delta_threshold}

  test_loom = scanpy.read_loom(
    test_loom_file, obs_names="cell_names", var_names="gene_names"
  )
  truth_loom = scanpy.read_loom(
    truth_loom_file, obs_names="cell_names", var_names="gene_names"
  )

  truth_cells = np.array(test_loom.X.sum(axis=1)).flatten()
  test_cells = np.array(truth_loom.X.sum(axis=1)).flatten()

  differences = [
    1 for (truth, test) in zip(truth_cells, test_cells) if (truth - test) != 0
  ]

  delta = len(differences) / len(truth_cells)

  if delta < threshold:
      sys.stdout.write(
          f"Matrices are identical: delta: {delta} delta_cutoff: {threshold}"
      )
      sys.exit(0)
  else:
      sys.stderr.write(
          f"Matrices are NOT identical: delta: {delta} delta_cutoff: {threshold}"
      )
      sys.exit(1)

  CODE
  >>>

  runtime {
    docker: docker
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    preemptible: 3
  }
}

task CompareH5adFilesATAC {
  input {
    File truth_h5ad
    File test_h5ad
    String docker = "python:3.10.0-buster"
    Int disk_size_gb = ceil(size(truth_h5ad, "GiB") + size(test_h5ad, "GiB")) + 200
    Int memory_gb = 32
  }

  command <<<

    set -eo pipefail

    pip3 install anndata
    
    python3 <<CODE
    
    import anndata as ad
    import numpy as np
    import pandas as pd
    
    def compare_atac(test,truth):
        print(truth.obs)
        print(test.obs)
        truth.obs.describe()
        test.obs.describe()
        # Find the intersection of barcodes
        shared_indices = truth.obs.index.intersection(test.obs.index)
        truth_shared = truth[shared_indices]
        test_shared = test[shared_indices]
        # Look at length of barcodes and barcodes shared
        number_truth_barcodes=len(truth)
        print("Number of Truth barcodes: ", number_truth_barcodes)
        number_test_barcodes=len(test)
        print("Number of Test barcodes: ", number_test_barcodes)
        percent_barcodes_shared_truth = len(truth_shared)/len(truth)*100
        print("% Truth barcodes shared ", str(percent_barcodes_shared_truth))
        percent_barcodes_shared_test = len(test_shared)/len(test)*100
        print("% Test barcodes shared ", str(percent_barcodes_shared_test))
        stats=np.corrcoef(truth_shared.obs.n_fragment, y=test_shared.obs.n_fragment)
        print(stats)
        return stats
    
    truth_h5ad = "~{truth_h5ad}"
    test_h5ad = "~{test_h5ad}"
    truth = ad.read_h5ad(truth_h5ad)
    test = ad.read_h5ad(test_h5ad)
    
    truth_obs = pd.DataFrame(truth.obs)
    test_obs = pd.DataFrame(test.obs)
    
    print("Now running obs equivalence check")
    
    if truth_obs.equals(test_obs)==True:
        print("pass")
    else:
        print("Files are not identical, running additional checks")
        stats = compare_atac(test,truth)
        value=stats[0,1]
        if value>0.990:
            print("pass")
        else:
            exit("Files are not similar enough to pass test")
    
    print("Done running matrix equivalence check")
    
    CODE 
  >>>

  runtime {
    docker: docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_gb} GiB"
    preemptible: 3
  }
}

task CompareH5adFilesGEX {
  input {
    File truth_h5ad
    File test_h5ad
    String docker = "python:3.10.0-buster"
    Int disk_size_gb = ceil(size(truth_h5ad, "GiB") + size(test_h5ad, "GiB")) + 200
    Int memory_gb = 32
  }

  command <<<

    set -eo pipefail

    pip3 install anndata
    
    python3 <<CODE
    
    import anndata as ad
    import numpy as np
    import pandas as pd
    
    truth_h5ad = "~{truth_h5ad}"
    test_h5ad = "~{test_h5ad}"
    truth = ad.read_h5ad(truth_h5ad)
    test = ad.read_h5ad(test_h5ad)
    
    for x in truth.obs.columns:
        z = test.obs[x]
        y = truth.obs[x]
        if z.equals(y)==False:
            print("Cell Metric Column does not match:")
            print(x)
            print("Sum of test: ")
            print(z.sum())
            print("Sum of truth: ")
            print(y.sum())
            if x == "doublet_score":
                print("Doublet score is allowed to be different")
            else: 
                exit("Cell Metric does not match")
    print("Comparing test gene metrics to truth gene metrics using truth as ref")
    for x in truth.var.columns:
        z = test.var[x]
        y = truth.var[x]
        if z.equals(y)==False:
            print("Gene Metric Column does not match:")
            print(x)
    print("Making gene_names unique")
    test.var_names_make_unique()
    truth.var_names_make_unique()
    genes_correct=True
    for x in truth.var.columns:
        z = test.var[x]
        y = truth.var[x]
        if z.equals(y)==False:
            print("Gene metric does not match after making gene names unique")
            print(x)
            genes_correct=False
    print("Done")
    print("If no warning above Done, gene metrics match now that they are unique")

    print("Testing for new obs columns in test data set:")
    for x in test.obs.columns:
        if x not in truth.obs.columns:
            print("Column not in truth", x)
    print("Done")
    print("If no warning above Done, no new obs columns in test matrix")

    print("Testing for new var columns in test data set:")
    for x in test.var.columns:
        if x not in truth.var.columns:
            print("Column not in truth", x)
    print("Done")
    print("If no warning above Done, no new var columns in test matrix")
    print("Testing matrix count sums")
    if test.X.sum()==truth.X.sum():
        print("Counts match")
    else:
        print("Counts do not match")
        exit("Counts do not match")
    if genes_correct==False:
        exit("Gene metrics do not match")

    print("Done with equivalence check")
    
    CODE 
  >>>

  runtime {
    docker: docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_gb} GiB"
    preemptible: 3
  }
}

task CompareSnapTextFiles {
  input {
    Array[File] test_text_files
    Array[File] truth_text_files
  }

  command <<<
    exit_code=0

    test_files_length=~{length(test_text_files)}
    truth_files_length=~{length(truth_text_files)}

    while read -r a && read -r b <&3; do
      sort_files() {
       sort -t ',' -k2,2 -k3,3n -k4,4n "$1" | cut -d',' -f2,3,4 > "${1%.csv}.sorted.csv"
      }

      calc_md5() {
        md5sum "$1" | cut -d ' ' -f1
      }

       if [[ "$a" == *_fragments.csv && "$b" == *_fragments.csv ]] || [[ "$a" == *_binCounts_10000.csv && "$b" == *_binCounts_10000.csv ]]; then
         echo "Sorting File $a and $b"
         sort_files "$a"
         sort_files "$b"

         echo "Calculating md5sums for $a and $b"
         md5_a=$(calc_md5 "${a%.csv}.sorted.csv")
         md5_b=$(calc_md5 "${b%.csv}.sorted.csv")

         if [[ "$md5_a" == "$md5_b" ]]; then
             echo "Files $a and $b are identical"
         else
             echo "Files $a and $b are NOT identical"
             diff ${a%.csv}.sorted.csv ${b%.csv}.sorted.csv > diffs.txt
             exit_code=1
             echo "Diff between ${a%.csv}.sorted.csv  and ${b%.csv}.sorted.csv:" >&2
             cat diffs.txt >&2
         fi
       else
         echo "Sorting File $a and $b"
         sort "$a" > "${a%.csv}.sorted.csv"
         sort "$b" > "${b%.csv}.sorted.csv"

         echo "Calculating md5sums for ${a%.csv}.sorted.csv and ${b%.csv}.sorted.csv"
         md5_a=$(calc_md5 "${a%.csv}.sorted.csv")
         md5_b=$(calc_md5 "${b%.csv}.sorted.csv")

         if [ $md5_a = $md5_b ]; then
           echo "Files $a and $b are identical"
         else
           echo "Files ${a%.csv}.sorted.csv and ${b%.csv}.sorted.csv have different md5sums."
           diff ${a%.csv}.sorted.csv ${b%.csv}.sorted.csv > diffs.txt
           exit_code=1
           echo "Diff between ${a%.csv}.sorted.csv and ${b%.csv}.sorted.csv:" >&2
           cat diffs.txt >&2
         fi
       fi
    done < ~{write_lines(test_text_files)} 3<~{write_lines(truth_text_files)}

    exit $exit_code
  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 50 HDD"
    memory: "25 GiB"
    preemptible: 3
  }
}


task CompareLibraryFiles {
  input {
    File test_text_file
    File truth_text_file
  }

  command {
    exit_code=0

    a=~{test_text_file}
    b=~{truth_text_file}

    echo "Sorting files $a and $b"
    sort "$a" > "a.sorted"
    sort "$b" > "b.sorted"

    echo "Calculating md5sums for $a and $b"
    md5_a=$(md5sum "a.sorted" | cut -d ' ' -f1)
    md5_b=$(md5sum "b.sorted" | cut -d ' ' -f1)

    if [ $md5_a = $md5_b ]; then
      echo "Files $a.sorted and $b.sorted have matching md5sums and are the same."
    else
      echo "Files $a.sorted and $b.sorted have different md5sums."

      # Compare the files, excluding specific lines
      excluded_lines="percent_doublets|keeper_cells|keeper_mean_reads_per_cell|keeper_median_genes|percent_keeper|percent_usable"
        
      # Store the diff result, but only check non-excluded lines
      diff_output=$(diff <(grep -v -E $excluded_lines a.sorted) <(grep -v -E $excluded_lines b.sorted))

      if [ -z "$diff_output" ]; then
          echo "Files a.sorted and $b.sorted are the same when excluding specified lines."
      else
        echo "Files a.sorted and b.sorted have differences in non-excluded lines."
        echo "$diff_output"
        exit_code=2
      fi
    fi
    echo "Exiting with code $exit_code"
    exit $exit_code
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 100 HDD"
    memory: "50 GiB"
    preemptible: 3
  }
}

task CompareH5Files {
  input {
    File test_h5
    File truth_h5
  }

  command {

    apt update
    apt install -y hdf5-tools

    h5diff -v ~{test_h5} ~{truth_h5} > diff_output.txt
    diff_exit_code=$?

    echo "H5diff output:"
    # Print the diff output to the console
    cat diff_output.txt

    if [ $diff_exit_code -eq 0 ]; then
      echo "H5 files are identical."
    elif [ $diff_exit_code -eq 1 ]; then
      echo "H5 files differ."
    else
      echo "Error comparing H5 files."
    fi

    echo "Exiting with code $diff_exit_code"
    exit $diff_exit_code

  }

  runtime {
    docker: "ubuntu:20.04"
    disks: "local-disk 150 HDD"
    memory: "80 GiB"
    preemptible: 3
  }
}



