version 1.0

task ValidateSmartSeq2Plate {
  input {
    File? loom_output
    File truth_loom

    Int disk_size = ceil(size(loom_output,"GiB") + size(truth_loom, "GiB") + 10)
  }

  command <<<

    # catch intermittent failures
    set -eo pipefail

   python3 /tools/loomCompare.py --truth-loom ~{truth_loom} --check-loom ~{loom_output} --delta-cutoff 10

  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-loom-output:0.0.3-fk-2"
    cpu: 1
    memory: "8 GiB"
    disks: "local-disk 1${disk_size} HDD"
  }
}

task CompareMetrics {
    input {
      File? target_metrics
      String expected_metrics_hash
    }
    Boolean target_matrics_defined = if defined(target_metrics) then true else false

  command <<<
    # catch intermittent failures
    set -eo pipefail

   if ~{target_matrics_defined}; then
     # this parses the picard metrics file with awk to remove all the run-specific comment lines (#)
     target_metrics_hash=$(cat "~{target_metrics}" | awk 'NF && $1!~/^#/' | md5sum | awk '{print $1}')

     if [ "$target_metrics_hash" != "~{expected_metrics_hash}" ]; then
       >&2 echo "target_metrics_hash ($target_metrics_hash) did not match expected hash (${expected_metrics_hash})"
       fail=true
     fi
   else
     >&2 echo "target_metrics_hash is not available"
   fi

   if [ $fail == "true" ]; then exit 1; fi

    >>>
  runtime {
    docker: "ubuntu:16.04"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
  }
}

task CompareCounts {
    input {
      String counts
      String expected_counts_hash
    }

  command <<<

    # catch intermittent failures
    set -eo pipefail

    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed. We parse the first 7 columns because a bug in RSEM
    # makes later columns non-deterministic.
    #counts_hash=$(awk 'NR>2' "~{counts}" | md5sum | awk '{print $1}')
    counts_hash=$counts
    if [ "$counts" != "$expected_counts_hash" ]; then
        echo "Strings are not equal"
        fail=true
    fi

    if [ "$fail" == "true" ]; then exit 1; fi
  >>>

  runtime {
    docker: "ubuntu:16.04"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
  }
}


