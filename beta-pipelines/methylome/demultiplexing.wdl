version 1.0

workflow Demultiplexing {

  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id
    File genome_fa
    File genome_ver
  }


  scatter(idx in range(length(fastq_input_read1))) {
    call Demultiplex {
      input:
        fastq_input_read1 = [fastq_input_read1[idx]],
        fastq_input_read2 = [fastq_input_read2[idx]],
        random_primer_indexes = random_primer_indexes,
        plate_id = plate_id
    }
  }
 call hisat_build_methyl_index {
   input:
     genome_fa = genome_fa,
     genome_ver = genome_ver
     }
  call Mapping {
    input:
      tarred_demultiplexed_fastqs = Demultiplex.output_fastqs,
      genome_fa = genome_fa,
      genome_ver = genome_ver,
      index_files = hisat_build_methyl_index.output_index_files
  }
  output {
    Array[File] output_fastqs = Demultiplex.output_fastqs
  }
}

task Demultiplex {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id

    String docker_image = "ekiernan/yap_hisat:v4"
    Int disk_size = 50
    Int mem_size = 10
  }

  command <<<
    set -euo pipefail
    declare -a fastq1_file=(~{sep=' ' fastq_input_read1})
    declare -a fastq2_file=(~{sep=' ' fastq_input_read2})

    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels \
    -g file:~{random_primer_indexes} \
    -o ~{plate_id}-{name}-R1.fq.gz \
    -p ~{plate_id}-{name}-R2.fq.gz \
    $fastq1_file \
    $fastq2_file \
    > ~{plate_id}.stats.txt

    # remove the fastq files that end in unknown-R1.fq.gz and unknown-R2.fq.gz
    rm *-unknown-R{1,2}.fq.gz

    # zip up all the output fq.gz files
    tar -zcvf ~{plate_id}.cutadapt_output_files.tar.gz *.fq.gz
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    File output_fastqs = "~{plate_id}.cutadapt_output_files.tar.gz"
    File stats = "~{plate_id}.stats.txt"
  }
}


task hisat_build_methyl_index {
  input {
  File genome_fa
  String genome_ver
  Int memory = 50
  Int disk = 1000
 # String hisat_command = "hisat-3n-build --base-change C,T --repeat-index"
}
command <<<
  set -euo pipefail

  /hisat-3n/hisat-3n-build --base-change C,T --repeat-index ~{genome_fa} ~{genome_ver}

>>>
output {
  Array[File] output_index_files = glob("hg38*")
}
runtime {
  docker: "ekiernan/hisat3n-python:v1"
  memory: memory + " GiB"
  disks: "local-disk ~{disk} HDD"
  disk: disk + " GB" # TES
}
}

task Mapping {
  input {
    Array[File] tarred_demultiplexed_fastqs
    Array[File] index_files

    File genome_fa
    String genome_ver
    #String hisat_command = "hisat-3n-build --base-change C,T --repeat-index"

    String docker_image = "ekiernan/hisat3n-python:v1"
    Int disk_size = 50
    Int mem_size = 500
  }

  command <<<
    set -euo pipefail

    # /hisat-3n/hisat-3n-build --base-change C,T --repeat-index ~{genome_fa} ~{genome_ver}

    declare -a fastq_files=(~{sep=' ' tarred_demultiplexed_fastqs})

    # untar the demultiplexed fastq files
    tar -zxvf ${fastq_files}

    # put all R1 and R2 files into arrays
    r1=(*R1.fq.gz)
    r2=(*R2.fq.gz)


    for (( i=0; i<${#r1[@]}; i++ )); do
    # define the output file names based on the input file names
    output_prefix=${r1[$i]%%.*} # remove the .r1.fq.gz suffix
    output_bam=${output_prefix}.bam

    /hisat-3n/hisat-3n \
    {config[hisat3n_dna_reference]} \
    -q \
    -1 ${r1[$i]} \
    -2 ${r2[$i]} \
    --directional-mapping-reverse \  # this can speed up 2X as the snmC reads are directional
    --base-change C,T \
    {repeat_index_flag} \
    --no-spliced-alignment \  # this is important for DNA mapping
    --no-temp-splicesite \
    -t \
    --new-summary \
    --summary-file ${output_prefix}.stats \
    --threads {threads} \
    | \
    samtools view \
    -b -q 0 -o ${output_bam}  # do not filter any reads in this step
    done

  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    Array[File] unsorted_bam_files = "output.bam"
  }
}
