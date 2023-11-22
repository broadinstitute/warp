version 1.0

task CompareVcfs {
  input {
    File file1
    File file2
    String patternForLinesToExcludeFromComparison = ""
  }

  command {
    set -eo pipefail

    if [ -z ~{patternForLinesToExcludeFromComparison} ]; then
      diff <(gunzip -c -f ~{file1}) <(gunzip -c -f ~{file2})
    else
      echo "It's defined!"
      diff <(gunzip -c -f ~{file1} | grep -v '~{patternForLinesToExcludeFromComparison}') <(gunzip -c -f ~{file2} | grep -v '~{patternForLinesToExcludeFromComparison}')
    fi
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
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
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
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
    gatk VCFComparator -R ~{ref_fasta}  -V:actual ~{actual} -V:expected ~{expected} ~{extra_args} ~{if(warn_on_error) then "--warn-on-errors" else ""} --finish-before-failing
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/gatk-vcfcomparator@sha256:4c1b32dd89c46af52e68ae34f99db483ba07b08def2479d145a185de0b2d9a4a"
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
        diff $a.sorted $b.sorted > diffs.txt
        exit_code=1
        echo "Diff between $a.sorted and $b.sorted:" >&2
        cat diffs.txt >&2
      fi

      # catting the diffs.txt on STDOUT as that's what's expected.
      cat diffs.txt

    done < ~{write_lines(test_text_files)} 3<~{write_lines(truth_text_files)}

    exit $exit_code
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
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
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
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
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
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
  Int memory_mb = 500000
  Int java_memory_size = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command {
    set -e
    set -o pipefail

    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
    CompareSAMs \
          ~{test_bam} \
          ~{truth_bam} \
          O=comparison.tsv \
          LENIENT_HEADER=~{lenient_header} \
          LENIENT_LOW_MQ_ALIGNMENT=~{lenient_low_mq}
  }

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
    diff <(gunzip -c -f ~{test_zip}) <(gunzip -c -f ~{truth_zip})
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
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
        exit("Files are not identical")
    
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
    
    truth_obs = pd.DataFrame(truth.obs)
    test_obs = pd.DataFrame(test.obs)
    
    truth_var = pd.DataFrame(truth.var)
    test_var = pd.DataFrame(test.var)
    
    truth_sum = truth.X.sum()
    test_sum = test.X.sum()
    
    print("Now running equivalence check")
    
    if truth_obs.equals(test_obs)==True and truth_var.equals(test_var)==True and truth_sum==test_sum:
        print("pass")
    else:
        exit("Files are not identical")
    
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
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
    memory: "25 GiB"
    preemptible: 3
  }
}


