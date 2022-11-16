version 1.0

#WORKFLOW DEFINITION   
workflow VUMCCramQC {
    input {
        Array[File] input_crams
        File ref_fasta
    }
  
  String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
  String gatk_path = "/gatk/gatk"

    output {
      Array[File] validation_reports = ValidateCRAM.validation_report
    }

  scatter (oneInput in input_crams) {
    String sample_basename = basename(oneInput, ".cram")
    call ValidateCRAM {
        input:
            input_cram = oneInput,
            sample_basename = sample_basename,
            docker = gatk_docker,
            gatk_path = gatk_path
    }
  }
}

# TASK DEFINITIONS
# Validate a cram using Picard ValidateSamFile
task ValidateCRAM {
  input {
    # Command parameters
    File input_cram
    File ref_fasta
    String sample_basename
    String? validation_mode
    String gatk_path
  
    # Runtime parameters
    String docker
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50
  }
    
  Int disk_size = ceil(size(input_cram, "GB")) + addtional_disk_space_gb
  String output_name = "${sample_basename}_${validation_mode}.txt"
 
  command {
    ${gatk_path} \
      ValidateSamFile \
      --INPUT ${input_cram} \
      --OUTPUT ${output_name} \
      --MODE ${default="VERBOSE" validation_mode}
      --REFERENCE_SEQUENCE ${ref_fasta}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File validation_report = "${output_name}"
  }
}