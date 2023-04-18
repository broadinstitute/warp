version 1.0

#WORKFLOW DEFINITION   
workflow VUMCCramQC {
  input {
    Array[File] input_crams
    String sample_name
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    String gatk_path = "/gatk/gatk"
    String samtools_docker = "staphb/samtools:latest"
    File? reference_file
    File? reference_file_dict
    File? reference_file_fai
  }

  scatter (input_cram in input_crams) {
    call ValidateCRAM {
      input:
        input_cram = input_cram,
        docker = gatk_docker,
        gatk_path = gatk_path,
        reference_file = reference_file,
        reference_file_dict = reference_file_dict,
        reference_file_fai = reference_file_fai
    }
  }

  call SummerizeQC {
    input:
      qc_reports = ValidateCRAM.validation_report,
      qc_codes = ValidateCRAM.cram_qc_code,
      sample_name = sample_name,
  }

  if(SummerizeQC.qc_failed==0) {
    scatter (input_cram in input_crams){
      call CountCRAM {
        input:
          input_cram = input_cram,
          docker = samtools_docker,
      }
    }

    call SumUp {
      input: 
        sample_name = sample_name,
        mapped_files = CountCRAM.mapped_file,
        unmapped_files = CountCRAM.unmapped_file,
    }
  }

  output {
    File cram_qc_reports = SummerizeQC.qc_reports
    Int cram_qc_failed = SummerizeQC.qc_failed
    Int? unmapped_reads = SumUp.NumberUnmappedReads
    Int? mapped_reads = SumUp.NumberMappedReads
  }
}

# TASK DEFINITIONS
# Validate a cram using Picard ValidateSamFile
task ValidateCRAM {
  input {
    File input_cram
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

  Int disk_size = ceil(size(input_cram, "GB")) + addtional_disk_space_gb
  String output_name = "summary.txt"
  String res_file = "res.txt"

  command <<<
    ~{gatk_path} \
      ValidateSamFile \
      --INPUT ~{input_cram} \
      --OUTPUT validate.summary  ~{"--REFERENCE_SEQUENCE " + reference_file} \
      --MODE ~{validation_mode} \
      --IGNORE MISSING_TAG_NM --IGNORE MATE_NOT_FOUND

    status=$?
    echo "$status" > ~{res_file}

    f="$(basename -- ~{input_cram})"
    echo "$f" > ~{output_name}
    if [[ -s validate.summary ]]; then
      cat validate.summary >> ~{output_name}
      rm -f validate.summary
    else
      echo "no summary genereated" >> ~{output_name}
    fi
  >>>

  runtime {
    docker: docker
    preemptible: 3
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File validation_report = "~{output_name}"
    Int cram_qc_code = read_int("~{res_file}")
  }
}

task SummerizeQC {
  input {
    Array[File] qc_reports
    Array[Int] qc_codes
    String sample_name
  }
    
  String output_name = "~{sample_name}.txt"
  String res_file = "res.txt"
  String toRunCount = "logical_runif.txt"

  command <<<
    cat ~{sep=" " qc_reports} > ~{output_name}
    python3 -c 'print(~{sep="+" qc_codes})' > ~{res_file}

  >>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: 3
    memory: "4 GB"
    disks: "local-disk 10 HDD"
  }

  output {
    File qc_reports = "~{output_name}"
    Int qc_failed = read_int("~{res_file}")
  }
}

task CountCRAM{
  input{
    # Command parameters
    File? input_cram

    # Runtime parameters
    String docker
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50
  }

  Int disk_size = ceil(size(input_cram, "GB")) + addtional_disk_space_gb
  String NumUnmapped = "Unmapped.txt"
  String NumMapped = "Mapped.txt"

  command <<<
    samtools flagstat ~{input_cram} | grep "mapped (" | grep -v "primary" | cut -f1 -d' ' > ~{NumMapped}

    samtools view -c -f4 ~{input_cram} > ~{NumUnmapped}
  >>>

  runtime{
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output{
    File unmapped_file = "~{NumUnmapped}"
    File mapped_file = "~{NumMapped}"
  }
}

task SumUp{
  input{
    # Command parameters
    String sample_name
    Array[File]? unmapped_files
    Array[File]? mapped_files
  }

  String FinalNumUnmapped = "${sample_name}_final_Unmapped.txt"
  String FinalNumMapped = "${sample_name}_final_Mapped.txt"

  command <<<
    # for array of file, we need to use sep to join them
    cat ~{sep=" " mapped_files} | awk '{ sum += $1 } END { print sum }' > ~{FinalNumMapped}
    cat ~{sep=" " unmapped_files} | awk '{ sum += $1 } END { print sum }' > ~{FinalNumUnmapped}
  >>>

  runtime{
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    memory: "2 GB"
    disks: "local-disk 5 HDD"
  }
  output{
    Int NumberUnmappedReads = read_int("~{FinalNumUnmapped}")
    Int NumberMappedReads = read_int("~{FinalNumMapped}")
  }
}
