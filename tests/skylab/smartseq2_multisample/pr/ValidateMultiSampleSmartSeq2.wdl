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

   python3 /usr/gitc/loomCompare.py --truth-loom ~{truth_loom} --check-loom ~{loom_output} --delta-cutoff 10

  >>>
  
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/pytools:1.0.0-1661263730"
    cpu: 1
    memory: "8 GiB"
    disks: "local-disk 1${disk_size} HDD"
  }
}
