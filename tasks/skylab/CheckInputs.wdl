version 1.0

task checkInputArrays {
  input {
    Boolean paired_end
    Array[String] input_ids
    Array[String]? input_names
    Array[String] fastq1_input_files
    Array[String] fastq2_input_files
  }
  Int len_input_ids = length(input_ids)
  Int len_fastq1_input_files = length(fastq1_input_files)
  Int len_fastq2_input_files = length(fastq2_input_files)
  Int len_input_names = if defined(input_names) then length(select_first([input_names])) else 0

  meta {
    description: "checks input arrays to ensure that all arrays are the same length"
  }

  command {
    set -e

    if [[ ~{len_input_ids} != ~{len_fastq1_input_files} ]]
      then
      echo "ERROR: Different number of arguments for input_id and fastq1 files"
      exit 1;
    fi

    if [[ ~{len_input_names} != 0 && ~{len_input_ids} != ~{len_input_names} ]]
        then
        echo "ERROR: Different number of arguments for input_name and input_id"
        exit 1;
    fi

    if  ~{paired_end} && [[ ~{len_fastq2_input_files} != ~{len_input_ids} ]]
      then
      echo "ERROR: Different number of arguments for sample names and fastq1 files"
      exit 1;
    fi
    exit 0;
  }

  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "1 GiB"
    disks: "local-disk 1 HDD"
  }

}
