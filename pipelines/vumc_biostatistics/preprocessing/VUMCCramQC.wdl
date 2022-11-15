version 1.0
#Workflow definition   
workflow VUMCCramQC {
    input {
        File input_cram
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        String output_basename
        String sample_names

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
        String gatk_path = "/gatk/gatk"
        String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
        String samtools_path = "samtools"
    }

    output {
    Array[File] validation_reports = ValidateBam.validation_report
    }

    #Need to make into a bam file since cram is incompatible with ValidateSamFile
    String sample_basename = basename(input_cram, ".cram")
    call CramToBamTask {
      input:
        input_cram = input_cram,
        sample_name = sample_basename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = gitc_docker,
        samtools_path = samtools_path
    }

    call ValidateSamFile {
        input:
            input_bam = CramToBamTask.output_bam
            docker = gatk_docker,
            gatk_path = gatk_path
    }
}

# TASK DEFINITIONS
task CramToBamTask {
  input {
    # Command parameters
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_cram
    String sample_name

    # Runtime parameters
    String docker
    Int? machine_mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
    String samtools_path
  }
    Float output_bam_size = size(input_cram, "GB") / 0.60
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 20
  
  command {
    set -e
    set -o pipefail

    ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram} |
    ~{samtools_path} view -b -o ~{sample_name}.bam -
    ~{samtools_path} index -b ~{sample_name}.bam
    mv ~{sample_name}.bam.bai ~{sample_name}.bai
  }
  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 15]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
 }
  output {
    File output_bam = "~{sample_name}.bam"
    File output_bai = "~{sample_name}.bai"
  }
}

# Validate a SAM or BAM using Picard ValidateSamFile
task ValidateBAM {
  input {
    # Command parameters
    File input_bam
    String sample_basename
    String? validation_mode
    String gatk_path
  
    # Runtime parameters
    String docker
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50
  }
    
  Int disk_size = ceil(size(input_bam, "GB")) + addtional_disk_space_gb
  String output_name = "${sample_basename}_${validation_mode}.txt"
 
  command {
    ${gatk_path} \
      ValidateSamFile \
      --INPUT ${input_bam} \
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

