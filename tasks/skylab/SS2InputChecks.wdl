version 1.0

task checkSS2Input {
  input {
    File fastq1
    File? fastq2
    Boolean paired_end
    Boolean force_no_check

    Int disk = ceil((size(fastq1, "GiB") + size(fastq2, "GiB")) + 10)
    String paired_end_str = if paired_end then "true" else "false"
  }

  meta {
    description: "checks ss2 input values and fails the pipeline immediately"
  }

  command {
    set -e

    ## Set pass to true
    pass="true"

    ## Perform checks
    if [[ ! -f "${fastq1}" ]]
    then
    pass="false"
    echo "ERROR: fastq1 (${fastq1}) not provided"
    fi

    if [[ ${paired_end_str} == "true" && ! -f "${fastq2}" ]]
    then
    pass="false"
    echo "ERROR: paired_end is set and fastq2 is not provided"
    fi

    if [[ ${force_no_check} == "true" ]]
    then
       echo "force_no_check is set: Ignoring input checks"
       exit 0;
    fi

    ## fail if any tests failed, ignore if force_no_check is set
    if [[ $pass == "true" ]]
    then
      exit 0;
    else
      exit 1;
    fi
  }

  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "1 GiB"
    disks: "local-disk ${disk} HDD"
  }

}
