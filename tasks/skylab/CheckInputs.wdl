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
  Int disk = 10

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
    docker: "bashell/alpine-bash@sha256:965a718a07c700a5204c77e391961edee37477634ce2f9cf652a8e4c2db858ff"
    cpu: 1
    memory: "1 GiB"
    disks: "local-disk ~{disk} HDD"
    disk: disk + " GB" # TES
  }

}


task checkOptimusInput {
  input {
    String cloud_provider
    #String SAS_TOKEN
    File r1_fastq
    String counting_mode
    Boolean force_no_check
    Boolean count_exons
    Int disk = ceil(size(r1_fastq, "GiB")) + 50
    Int machine_mem_mb = 1000
    Int cpu = 1
    Int tenx_chemistry_version
    String gcp_whitelist_v2
    String gcp_whitelist_v3
    String azure_whitelist_v2
    String azure_whitelist_v3
    Boolean ignore_r1_read_length
    String alpine_docker_path
  }  

  meta {
    description: "checks optimus input values and fails the pipeline immediately"
  }

  command <<<
    set -e

    ## Set pass to true
    pass="true"
    
    ## Need to gunzip the r1_fastq
    zcat ~{r1_fastq} | head -n2 > r1.fastq
    FASTQ=r1.fastq
    echo 'this is the fastq:' $FASTQ
    R1=$(awk 'NR==2' $FASTQ)
    COUNT=$(echo ${#R1})
    echo 'this is the read:' $R1
    echo 'this is the UMI/barcode count:' $COUNT

    ## Perform checks
    if [[ ! ("~{counting_mode}" == "sc_rna" || "~{counting_mode}" == "sn_rna") ]]
    then
      pass="false"
      echo "ERROR: Invalid value \"${counting_mode}\" for input \"counting_mode\""
    fi

    if [[ ~{force_no_check} == "true" ]]
    then
       echo "force_no_check is set: Ignoring input checks"
       exit 0;
    fi

    if [[ "~{counting_mode}" == "sc_rna" ]]
    then
      if [[ ~{count_exons} == "true" ]]
      then
        pass="false"
        echo "ERROR: Invalid value count_exons should not be used with \"${counting_mode}\" input."
      fi
    fi

    # Check for chemistry version to produce read structure and whitelist
    if [[ ~{tenx_chemistry_version} == 2 ]]
      then
      if [[ "~{cloud_provider}" == "gcp" ]]
      then
        WHITELIST="~{gcp_whitelist_v2}"
      elif [[ "~{cloud_provider}" == "azure" ]]
      then
        WHITELIST="~{azure_whitelist_v2}"
      else
        pass="false"
        echo "ERROR: Cloud provider must be either gcp or azure"
      fi
      echo "WHITELIST:" $WHITELIST
      echo $WHITELIST > whitelist.txt
      echo 16C10M > read_struct.txt
    elif [[ ~{tenx_chemistry_version} == 3 ]]
      then
      if [[ "~{cloud_provider}" == "gcp" ]]
      then
        WHITELIST="~{gcp_whitelist_v3}"
      elif [[ "~{cloud_provider}" == "azure" ]]
      then
        WHITELIST="~{azure_whitelist_v3}"
      else
        pass="false"
        echo "ERROR: Cloud provider must be either gcp or azure"
      fi
      echo "WHITELIST:" $WHITELIST
      echo $WHITELIST > whitelist.txt
      echo 16C12M > read_struct.txt
    else
      pass="false"
      echo "ERROR: Chemistry version must be either 2 or 3"
    fi
    
    if [[ ~{tenx_chemistry_version} == 2 && $COUNT != 26 && ~{ignore_r1_read_length} == "false" ]]
      then
      pass="false"
      echo "Read1 FASTQ does not match v2 chemistry; to override set ignore_r1_read_length to true"
    elif [[ ~{tenx_chemistry_version} == 3 && $COUNT != 28 && ~{ignore_r1_read_length} == "false" ]]
      then
      pass="false"
      echo "Read1 FASTQ does not match v3 chemistry; to override set ignore_r1_read_length to true"
    else
      pass="true"
    fi


    ## fail if any tests failed, ignore if force_no_check is set
    if [[ $pass == "true" ]]
    then
      exit 0;
    else
      exit 1;
    fi

    exit 0;
  >>>

  output {
    String whitelist_out = read_string("whitelist.txt")
    String read_struct_out = read_string("read_struct.txt")
  }
  runtime {
    docker: alpine_docker_path
    cpu: cpu
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    disk: disk + " GB" # TES
  } 
}
