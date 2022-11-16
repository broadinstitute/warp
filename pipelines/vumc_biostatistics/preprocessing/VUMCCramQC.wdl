version 1.0

#WORKFLOW DEFINITION   
workflow VUMCCramQC {
    input {
        Array[File] input_cram
        String sample_name
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
        String gatk_path = "/gatk/gatk"
    }

    output {
      File validation_reports = ValidateCRAM.validation_report
    }

    call ValidateCRAM {
        input:
            input_cram = input_cram,
            sample_basename = sample_name,
            docker = gatk_docker,
            gatk_path = gatk_path
    }
}

# TASK DEFINITIONS
# Validate a cram using Picard ValidateSamFile
task ValidateCRAM {
  input {
    # Command parameters
    Array[File] input_cram
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
      --MODE ${default="SUMMARY" validation_mode}
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