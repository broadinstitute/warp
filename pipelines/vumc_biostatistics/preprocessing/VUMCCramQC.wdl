version 1.0

#WORKFLOW DEFINITION   
workflow VUMCCramQC {
  input {
    Array[File] input_crams
    Array[File]? input_crais
    String sample_name
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    String gatk_path = "/gatk/gatk"
    File? reference_file
    File? reference_file_dict
    File? reference_file_fai
  }

  output {
    File cram_qc_reports = ValidateCRAM.validation_report
    Int cram_qc_failed = ValidateCRAM.cram_qc_failed
  }

  call ValidateCRAM {
    input:
      input_crams = input_crams,
      input_crais = input_crais,
      sample_name = sample_name,
      docker = gatk_docker,
      gatk_path = gatk_path,
      reference_file = reference_file,
      reference_file_dict = reference_file_dict,
      reference_file_fai = reference_file_fai
  }
}

# TASK DEFINITIONS
# Validate a cram using Picard ValidateSamFile
task ValidateCRAM {
  input {
    # Command parameters
    Array[File] input_crams
    Array[File]? input_crais
    String sample_name
    String validation_mode = "SUMMARY"
    String gatk_path
    File? reference_file
    File? reference_file_dict
    File? reference_file_fai
  
    # Runtime parameters
    String docker
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50
  }
    
  Int disk_size = ceil(size(input_crams, "GB")) + addtional_disk_space_gb
  String output_name = "${sample_name}_${validation_mode}.txt"
  String res_file = "${sample_name}_res.txt"
 
  command <<<
    echo "0" > ~{res_file}

    for input_cram in ~{sep=" " input_crams}
    do
      f="$(basename -- $input_cram)"

      ~{gatk_path} \
        ValidateSamFile \
        --INPUT $input_cram \
        --OUTPUT validate.summary  ~{"--REFERENCE_SEQUENCE " + reference_file} \
        --MODE ~{validation_mode} \
        --IGNORE MISSING_TAG_NM --IGNORE MATE_NOT_FOUND

      status=$?
      if [[ $status != 0 ]]; then
        echo "$status" > ~{res_file}
      fi

      echo "$f" >> ~{output_name}
      if [[ -s validate.summary ]]; then
        cat validate.summary >> ~{output_name}
        rm -f validate.summary
      else
        echo "no summary genereated" >> ~{output_name}
      fi
    done
  >>>

  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File validation_report = "~{output_name}"
    Int cram_qc_failed = read_int("~{res_file}")
  }
}