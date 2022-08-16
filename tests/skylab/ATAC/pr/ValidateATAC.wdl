version 1.0

task ValidateATAC {
  input {
      File snap
      File snapqc
      File bam

      Int required_disk = ceil((size(bam, "G") + size(snap, "G")) * 1.2)

      String expected_snap_hash
      String expected_snapqc_hash
      String expected_bam_hash
   }

  command <<<
    # catch intermittent failures
    set -eo pipefail

    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed. gzipped files are unzipped to avoid hashing malleable
    # metadata

    snap_hash=$(cat ~{snap} | md5sum |  awk '{print $1}')
    snapqc_hash=$(cat ~{snapqc} | md5sum |  awk '{print $1}')

    # calculate hash as above, but ignore run-specific bam headers
    bam_hash=$(samtools view ~{bam} | md5sum | awk '{print $1}')

    # test each output for equality, echoing any failure states to stdout
    fail=false
    if [ "$bam_hash" != "~{expected_bam_hash}" ]; then
      >&2 echo "bam_hash ($bam_hash) did not match expected hash (~{expected_bam_hash})"
      fail=true
    fi

    ## This is disabled because we need to use the breakout step results instead
    ## The snap file contains a time stamp in the header and therefore direct chksum 
    ## comparisons are not possible

    #if [ "$snap_hash" != "~{expected_snap_hash}" ]; then
    #  >&2 echo "snap_hash ($snap_hash) did not match expected hash (~{expected_snap_hash})"
    #  fail=true
    #fi

    if [ "$snapqc_hash" != "~{expected_snapqc_hash}" ]; then
      >&2 echo "snapqc_hash ($snapqc_hash) did not match expected hash (~{expected_snapqc_hash})"
      fail=true
    fi

    if [ $fail == "true" ]; then exit 1; fi
  >>>
  
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${required_disk} HDD"
  }
}
