version 1.0

task ValidateSmartSeq2SingleCell {
    input {
      File counts
      String expected_counts_hash

      File? target_metrics
      String expected_metrics_hash
    }
    Boolean target_matrics_defined = if defined(target_metrics) then true else false
    
  command <<<

    # catch intermittent failures
    set -eo pipefail

    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed. We parse the first 7 columns because a bug in RSEM
    # makes later columns non-deterministic.
    counts_hash=$(cut -f 1-7 "~{counts}" | md5sum | awk '{print $1}')

    if [ "$counts_hash" != "~{expected_counts_hash}" ]; then
      >&2 echo "counts_hash ($counts_hash) did not match expected hash (${expected_counts_hash})"
      fail=true
    fi

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
