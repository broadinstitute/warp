version 1.0

task HISAT2PairedEnd {
  input {
    File hisat2_ref
    File fastq1
    File? fastq2
    String ref_name
    String output_basename
    String sample_name

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
  Int machine_mem_mb = 16500
  Int cpu = 4
  # Using (fastq1 + fastq2) x 100 gives factor of a few buffer. BAM can be up to ~5 x (fastq1 + fastq2).
  # Need room for unsorted + sorted bam + temp sorting space + zipped and unzipped ref. Add 10 GiB buffer.
  Int disk = ceil((size(fastq1, "GiB") + size(fastq2, "GiB")) * 100 + size(hisat2_ref, "GiB") * 2 + 10)
  Int preemptible = 5
}
  meta {
    description: "HISAT2 alignment task will align paired-end fastq reads to reference genome."
  }

  parameter_meta {
    hisat2_ref: "HISAT2 reference"
    fastq1: "gz forward fastq file"
    fastq2: "gz reverse fastq file"
    ref_name: "the basename of the index for the reference genome"
    output_basename: "basename used for output files"
    sample_name: "sample name of input"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    # Note that files MUST be gzipped or the module will not function properly
    # This will be addressed in the future either by a change in how Hisat2 functions or a more
    # robust test for compression type.

    set -e

    # fix names if necessary.
    if [ "${fastq1}" != *.fastq.gz ]; then
        FQ1=${fastq1}.fastq.gz
        mv ${fastq1} ${fastq1}.fastq.gz
    else
        FQ1=${fastq1}
    fi
    if [ "${fastq2}" != *.fastq.gz ]; then
        FQ2=${fastq2}.fastq.gz
        mv ${fastq2} ${fastq2}.fastq.gz
    else
        FQ2=${fastq2}
    fi

    tar --no-same-owner -xvf "${hisat2_ref}"

    # run HISAT2 to genome reference with dedault parameters
    # --seed to fix pseudo-random number and in order to produce deterministics results
    # --secondary reports secondary alignments for multimapping reads. -k 10
    # searches for up to 10 primary alignments for each read
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -1 $FQ1 \
      -2 $FQ2 \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file ${output_basename}.log \
      --met-file ${output_basename}.hisat2.met.txt --met 5 \
      --seed 12345 \
      -k 10 \
      --secondary \
      -p ${cpu} -S >(samtools view -1 -h -o ${output_basename}_unsorted.bam)
    samtools sort -@ ${cpu} -O bam -o "${output_basename}.bam" "${output_basename}_unsorted.bam"
    samtools index "${output_basename}.bam"
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file = "${output_basename}.log"
    File met_file = "${output_basename}.hisat2.met.txt"
    File output_bam = "${output_basename}.bam"
    File bam_index = "${output_basename}.bam.bai"
  }
}

task HISAT2RSEM {
  input {
    File hisat2_ref
    File fastq1
    File? fastq2
    String ref_name
    String output_basename
    String sample_name

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    Int machine_mem_mb = 16500
    Int cpu = 4
    # Using (fastq1 + fastq2) x 100 gives factor of a few buffer. BAM can be up to ~5 x (fastq1 + fastq2).
    # Need room for unsorted + sorted bam + temp sorting space + zipped and unzipped ref. Add 10 GiB buffer.
    Int disk = ceil((size(fastq1, "GiB") + size(fastq2, "GiB")) * 100 + size(hisat2_ref, "GiB") * 2 + 10)
    Int preemptible = 5
  }

  meta {
    description: "This HISAT2 alignment task will align paired-end fastq reads to transcriptome only. "
  }

  parameter_meta {
    hisat2_ref: "HISAT2 reference"
    fastq1: "gz forward fastq file"
    fastq2: "gz reverse fastq file"
    ref_name: "the basename of the index for the reference genome"
    output_basename: "basename used for output files"
    sample_name: "sample name of input"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    # fix names if necessary.
    if [ "${fastq1}" != *.fastq.gz ]; then
        FQ1=${fastq1}.fastq.gz
        mv ${fastq1} ${fastq1}.fastq.gz
    else
        FQ1=${fastq1}
    fi

    if [ "${fastq2}" != *.fastq.gz ]; then
        FQ2=${fastq2}.fastq.gz
        mv ${fastq2} ${fastq2}.fastq.gz
    else
        FQ2=${fastq2}
    fi

    tar --no-same-owner -xvf "${hisat2_ref}"

    # increase gap alignment penalty to avoid gap alignment
    # --mp 1,1 --np 1 --score-min L,0,-0.1 are default parameters when rsem runs alignment by using bowtie2/Bowtie
    # --mp 1,1 and --np 1 will reduce mismatching penalty to 1 for all.
    # with no-splice-alignment no-softclip no-mixed options on, HISAT2 will only output concordant alignment without soft-cliping
    # --rdg 99999999,99999999 and --rfg 99999999,99999999 will set an infinite penalty to alignments with indels.
    # As a result, alignments with gaps or deletions are excluded.
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -1 $FQ1 \
      -2 $FQ2 \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file ${output_basename}.log \
      --met-file ${output_basename}.hisat2.met.txt --met 5 \
      -k 10 \
      --mp 1,1 \
      --np 1 \
      --score-min L,0,-0.1 \
      --secondary \
      --no-mixed \
      --no-softclip \
      --no-discordant \
      --rdg 99999999,99999999 \
      --rfg 99999999,99999999 \
      --no-spliced-alignment \
      --seed 12345 \
      -p ${cpu} -S >(samtools view -1 -h -o ${output_basename}.bam)
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file = "${output_basename}.log"
    File met_file = "${output_basename}.hisat2.met.txt"
    File output_bam = "${output_basename}.bam"
  }
}

task HISAT2SingleEnd {
input {
  File hisat2_ref
  File fastq
  String ref_name
  String output_basename
  String sample_name

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
  Int machine_mem_mb = 16500
  Int cpu = 4
  # Using fastq x 100 gives factor of a few buffer. BAM can be up to ~5 x fastq.
  # Need room for unsorted + sorted bam + temp sorting space + zipped and unzipped ref. Add 10 GiB buffer.
  Int disk = ceil((size(fastq, "GiB") * 100) + size(hisat2_ref, "GiB") * 2 + 10)
  Int preemptible = 5
}
  meta {
    description: "This HISAT2 alignment task will align single-end fastq reads to reference genome."
  }

  parameter_meta {
    hisat2_ref: "HISAT2 reference"
    fastq: "input fastq from single ended data"
    ref_name: "the basename of the index for the reference genome"
    output_basename: "basename used for output files"
    sample_name: "sample name of input"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e
    tar --no-same-owner -xvf "~{hisat2_ref}"

    # fix file names if necessary.
    if [ "~{fastq}" != *.fastq.gz ]; then
        FQ=~{fastq}.fastq.gz
        mv ~{fastq} ~{fastq}.fastq.gz
    else
        FQ=~{fastq}
    fi

    # The parameters for this task are copied from the HISAT2PairedEnd task.
    hisat2 -t \
      -x ~{ref_name}/~{ref_name} \
      -U $FQ \
      --rg-id=~{sample_name} --rg SM:~{sample_name} --rg LB:~{sample_name} \
      --rg PL:ILLUMINA --rg PU:~{sample_name} \
      --new-summary --summary-file "~{output_basename}.log" \
      --met-file ~{output_basename}.hisat2.met.txt --met 5 \
      --seed 12345 \
      -k 10 \
      --secondary \
      -p ~{cpu} -S >(samtools view -1 -h -o ~{output_basename}_unsorted.bam)
    samtools sort -@ ~{cpu} -O bam -o "~{output_basename}.bam" "~{output_basename}_unsorted.bam"
    samtools index "~{output_basename}.bam"
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file ="${output_basename}.log"
    File met_file ="${output_basename}.hisat2.met.txt"
    File output_bam = "${output_basename}.bam"
    File bam_index = "${output_basename}.bam.bai"
  }
}

task HISAT2InspectIndex {
  input {
    File hisat2_ref
    String ref_name

    # runtime values
    String docker =  "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    Int machine_mem_mb = 3850
    Int cpu = 1
    # use provided disk number or dynamically size on our own, with 10GiB of additional disk
    Int disk = ceil(size(hisat2_ref, "GiB") + 10)
    Int preemptible = 5
  }

  meta {
    description: "This task will test reference indexing files built for HISAT2 aligner."
  }

  parameter_meta {
    hisat2_ref: "HISAT2 reference"
    ref_name: "the basename of the index for the reference genome"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e
    tar --no-same-owner -xvf "${hisat2_ref}"
    hisat2-inspect --ss --snp \
       -s ${ref_name}/${ref_name} > hisat2_inspect.log
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file ="hisat2_inspect.log"
  }
}

task HISAT2RSEMSingleEnd {
  input {
    File hisat2_ref
    File fastq
    String ref_name
    String output_basename
    String sample_name

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    Int machine_mem_mb = 15000
    Int cpu = 4
    Int disk = ceil((size(fastq, "GiB")) * 100 + size(hisat2_ref, "GiB") * 2 + 10)
    Int preemptible = 5
  }

  meta {
    description: "This HISAT2 alignment task will align paired-end fastq reads to transcriptome only. "
  }

  parameter_meta {
    hisat2_ref: "HISAT2 reference"
    fastq: "gz fastq file"
    ref_name: "the basename of the index for the reference genome"
    output_basename: "basename used for output files"
    sample_name: "sample name of input"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    # fix names if necessary.
    if [ "${fastq}" != *.fastq.gz ]; then
        FQ=${fastq}.fastq.gz
        mv ${fastq} ${fastq}.fastq.gz
    else
        FQ=${fastq}
    fi


    tar --no-same-owner -xvf "${hisat2_ref}"

    # increase gap alignment penalty to avoid gap alignment
    # --mp 1,1 --np 1 --score-min L,0,-0.1 is default paramesters when rsem runs alignment by using bowtie2/Bowtie
    # --mp 1,1 and --np 1 will reduce mismatching penalty to 1 for all.
    # with no-splice-alignment no-softclip no-mixed options on, HISAT2 will only output concordant alignment without soft-cliping
    # --rdg 99999999,99999999 and --rfg 99999999,99999999 will set an infinite penalty to alignments with indels.
    # As a result, alignments with gaps or deletions are excluded.
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -U $FQ \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file ${output_basename}.log \
      --met-file ${output_basename}.hisat2.met.txt --met 5 \
      -k 10 \
      --mp 1,1 \
      --np 1 \
      --score-min L,0,-0.1 \
      --secondary \
      --no-mixed \
      --no-softclip \
      --no-discordant \
      --rdg 99999999,99999999 \
      --rfg 99999999,99999999 \
      --no-spliced-alignment \
      --seed 12345 \
      -p ${cpu} -S >(samtools view -1 -h -o ${output_basename}.bam)
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file = "${output_basename}.log"
    File met_file = "${output_basename}.hisat2.met.txt"
    File output_bam = "${output_basename}.bam"
  }
}
