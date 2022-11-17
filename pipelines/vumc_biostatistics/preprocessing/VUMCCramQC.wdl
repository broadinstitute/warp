version 1.0

#WORKFLOW DEFINITION   
workflow VUMCCramQC {
    input {
        File input_cram

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
        String gatk_path = "/gatk/gatk"
    }

    output {
      File cram_qc_reports = ValidateCRAM.validation_report
      Int cram_qc_failed = ValidateCRAM.cram_qc_failed
    }

    String sample_basename = basename(input_cram, ".cram")

    call ValidateCRAM {
        input:
            input_cram = input_cram,
            sample_basename = sample_basename,
            docker = gatk_docker,
            gatk_path = gatk_path
    }
}

# TASK DEFINITIONS
# Validate a cram using Picard ValidateSamFile
task ValidateCRAM {
  input {
    # Command parameters
    File input_cram
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
  String res_file = "${sample_basename}_res.txt"
 
  command {
    ${gatk_path} \
      ValidateSamFile \
      --INPUT ${input_cram} \
      --OUTPUT ${output_name} \
      --MODE ${default="SUMMARY" validation_mode}

    status=$?
    echo "$status" > ~{res_file}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File validation_report = "${output_name}"
    Int cram_qc_failed = read_int("~{res_file}")
  }
}