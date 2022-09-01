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
    docker: "bashell/alpine-bash:latest"
    cpu: 1
    memory: "1 GiB"
    disks: "local-disk 1 HDD"
  }

}


task checkOptimusInput {
  input {
    String chemistry
    String counting_mode
    Boolean force_no_check
    Boolean count_exons
    Int disk = 1
    Int machine_mem_mb = 1
    Int cpu = 1
  }  

  meta {
    description: "checks optimus input values and fails the pipeline immediately"
  }

  command {
    set -e

    ## Set pass to true
    pass="true"

    ## Perform checks
    if [[ ! ("${chemistry}" == "tenX_v2" || "${chemistry}" == "tenX_v3") ]]
    then
      pass="false"
      echo "ERROR: Invalid value \"${chemistry}\" for input \"chemistry\""
    fi

    if [[ ! ("${counting_mode}" == "sc_rna" || "${counting_mode}" == "sn_rna") ]]
    then
      pass="false"
      echo "ERROR: Invalid value \"${counting_mode}\" for input \"counting_mode\""
    fi

    if [[ ${force_no_check} == "true" ]]
    then
       echo "force_no_check is set: Ignoring input checks"
       exit 0;
    fi

    if [[ "${counting_mode}" == "sc_rna" ]]
    then
      if [[ ~{count_exons} == "true" ]]
      then
        pass="false"
        echo "ERROR: Invalid value count_exons should not be used with \"${counting_mode}\" input."
      fi
    fi

    ## fail if any tests failed, ignore if force_no_check is set
    if [[ $pass == "true" ]]
    then
      exit 0;
    else
      exit 1;
    fi

    exit 0;
  }

  runtime {
    docker: "bashell/alpine-bash:latest"
    cpu: cpu
    memory: "~{machine_mem_mb} GiB"
    disks: "local-disk ~{disk} HDD"
  }
  
}
