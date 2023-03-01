version 1.0

task StarAlignBamSingleEnd {
  input {
    File bam_input
    File tar_star_reference

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658781884"
    Int machine_mem_mb = ceil((size(tar_star_reference, "Gi")) + 6) * 1100
    Int cpu = 16
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 2.5) + (size(bam_input, "Gi") * 2.5))
    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 0
  }

  meta {
    description: "Aligns reads in bam_input to the reference genome in tar_star_reference"
  }

  parameter_meta {
    bam_input: "unaligned bam file containing genomic sequence, tagged with barcode information"
    tar_star_reference: "star reference tarball built against the species that the bam_input is derived from"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    # prepare reference
    mkdir genome_reference
    tar -xf "${tar_star_reference}" -C genome_reference --strip-components 1
    rm "${tar_star_reference}"

    STAR \
      --runMode alignReads \
      --runThreadN ${cpu} \
      --genomeDir genome_reference \
      --readFilesIn "${bam_input}" \
      --outSAMtype BAM Unsorted \
      --outSAMmultNmax -1 \
      --outSAMattributes All \
      --outSAMunmapped Within \
      --readFilesType SAM SE \
      --readFilesCommand samtools view -h \
      --runRNGseed 777
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} SSD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = "Aligned.out.bam"
    File alignment_log = "Log.final.out"
  }
}

task StarAlignFastqPairedEnd {
  input {
    File fastq1
    File fastq2
    File tar_star_reference

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658781884"
    Int machine_mem_mb = ceil((size(tar_star_reference, "Gi")) + 6) * 1100
    Int cpu = 16
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 2.5) + (size(fastq1, "Gi") * 5.0))
    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 3
  }

  meta {
    description: "Aligns reads in fastq1 and fastq2 to the reference genome in tar_star_reference"
  }

  parameter_meta {
    fastq1: "trimmed R1 FASTQ file containing genomic sequence"
    fastq2: "trimmed R2 FASTQ file containing genomic sequence"
    tar_star_reference: "star reference tarball built against the species that the bam_input is derived from"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    # prepare reference
    mkdir genome_reference
    tar -xf "${tar_star_reference}" -C genome_reference --strip-components 1
    rm "${tar_star_reference}"

    STAR \
    --genomeDir genome_reference \
    --runThreadN ${cpu} \
    --readFilesIn ~{fastq1} ~{fastq2} \
    --readFilesCommand "gunzip -c" \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --runRNGseed 777 \
    --limitBAMsortRAM 10000000000 \
    --quantMode GeneCounts
  }
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} SSD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File output_bam = "Aligned.sortedByCoord.out.bam"
  }

}

task StarAlignFastqMultisample {
  input {
    Array[File] fastq1_input_files
    Array[File] fastq2_input_files
    Array[String] input_ids
    File tar_star_reference

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658781884"
    Int machine_mem_mb = ceil((size(tar_star_reference, "Gi")) + 6) * 1100
    Int cpu = 16
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 2.5) + (size(fastq2_input_files, "Gi") * 2.0))
    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 3
  }

  meta {
    description: "Aligns reads in fastq1 and fastq2 to the reference genome in tar_star_reference"
  }

  parameter_meta {
    fastq1_input_files: "Array of trimmed R1 fastq files containing genomic sequence."
    fastq2_input_files: "Array of trimmed R2 fastq files containing genomic sequence."
    input_ids: "Array of input ids"
    tar_star_reference: "star reference tarball built against the species that the bam_input is derived from"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -e
    set -exo pipefail

    # prepare reference
    mkdir genome_reference
    tar -xf "~{tar_star_reference}" -C genome_reference --strip-components 1
    rm "~{tar_star_reference}"

    declare -a fastq1_files=(~{sep=' ' fastq1_input_files})
    declare -a fastq2_files=(~{sep=' ' fastq2_input_files})
    declare -a output_prefix=(~{sep=' ' input_ids})
    STAR --genomeLoad LoadAndExit --genomeDir genome_reference
    for (( i=0; i<${#output_prefix[@]}; ++i));
      do
        STAR \
          --genomeDir genome_reference \
          --runThreadN ~{cpu} \
          --readFilesIn ${fastq1_files[$i]} ${fastq2_files[$i]} \
          --readFilesCommand "gunzip -c" \
          --outSAMtype BAM SortedByCoordinate \
          --outReadsUnmapped Fastx \
          --runRNGseed 777 \
          --limitBAMsortRAM 10000000000 \
          --quantMode GeneCounts \
          --genomeLoad LoadAndKeep

        mv "Aligned.sortedByCoord.out.bam"   "${output_prefix[$i]}_Aligned.bam"
      done;
    STAR --genomeLoad Remove --genomeDir genome_reference
  >>>

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} SSD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    Array[File] output_bam = glob("*_Aligned.bam")
  }

}


task STARsoloFastq {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    File tar_star_reference
    File white_list
    Int chemistry
    String counting_mode
    String output_bam_basename
    Boolean? count_exons

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658781884"
    Int machine_mem_mb = 64000
    Int cpu = 8
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 3)) + ceil(size(r1_fastq, "Gi") * 20) +  ceil(size(r2_fastq, "Gi") * 20)
    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 3
  }

  meta {
    description: "Aligns reads in bam_input to the reference genome in tar_star_reference"
  }

  parameter_meta {
    r1_fastq: "input FASTQ file array"
    r2_fastq: "array of forward read FASTQ files"
    tar_star_reference: "star reference tarball built against the species that the bam_input is derived from"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -e

    UMILen=10
    CBLen=16
    if [ "~{chemistry}" == 2 ]
    then
        ## V2
        UMILen=10
        CBLen=16
    elif [ "~{chemistry}" == 3 ]
    then
        ## V3
        UMILen=12
        CBLen=16
    else
        echo Error: unknown chemistry value: "$chemistry". Should be one of "tenX_v2" or "texX_v3".
        exit 1;
    fi

    COUNTING_MODE=""
    if [[ "~{counting_mode}" == "sc_rna" ]]
    then
        ## single cell or whole cell
        COUNTING_MODE="Gene"
    elif [[ "~{counting_mode}" == "sn_rna" ]]
    then
    ## single nuclei
        COUNTING_MODE="GeneFull"
    else
        echo Error: unknown counting mode: "$counting_mode". Should be either sn_rna or sc_rna.
        exit 1;
    fi

    # prepare reference
    mkdir genome_reference
    tar -xf "~{tar_star_reference}" -C genome_reference --strip-components 1
    rm "~{tar_star_reference}"


    echo "UMI LEN " $UMILen
    if [[ ~{count_exons} ]]
    then
      STAR \
      --soloType Droplet \
      --soloStrand Unstranded \
      --runThreadN ~{cpu} \
      --genomeDir genome_reference \
      --readFilesIn "~{sep=',' r2_fastq}" "~{sep=',' r1_fastq}" \
      --readFilesCommand "gunzip -c" \
      --soloCBwhitelist ~{white_list} \
      --soloUMIlen $UMILen --soloCBlen $CBLen \
      --soloFeatures "Gene" \
      --clipAdapterType CellRanger4 \
      --outFilterScoreMin 30  \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloUMIdedup 1MM_Directional_UMItools \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes UB UR UY CR CB CY NH GX GN \
      --soloBarcodeReadLength 0
    fi

    STAR \
      --soloType Droplet \
      --soloStrand Unstranded \
      --runThreadN ~{cpu} \
      --genomeDir genome_reference \
      --readFilesIn "~{sep=',' r2_fastq}" "~{sep=',' r1_fastq}" \
      --readFilesCommand "gunzip -c" \
      --soloCBwhitelist ~{white_list} \
      --soloUMIlen $UMILen --soloCBlen $CBLen \
      --soloFeatures $COUNTING_MODE \
      --clipAdapterType CellRanger4 \
      --outFilterScoreMin 30  \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloUMIdedup 1MM_Directional_UMItools \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes UB UR UY CR CB CY NH GX GN \
      --soloBarcodeReadLength 0

    touch barcodes_sn_rna.tsv
    touch features_sn_rna.tsv
    touch matrix_sn_rna.mtx

    if [[ "~{counting_mode}" == "sc_rna" ]]
    then
      mv "Solo.out/Gene/raw/barcodes.tsv" barcodes.tsv
      mv "Solo.out/Gene/raw/features.tsv" features.tsv
      mv "Solo.out/Gene/raw/matrix.mtx"   matrix.mtx
    elif [[ "~{counting_mode}" == "sn_rna" ]]
    then
      if ! [[ ~{count_exons} ]]
      then
        mv "Solo.out/GeneFull/raw/barcodes.tsv" barcodes.tsv
        mv "Solo.out/GeneFull/raw/features.tsv" features.tsv
        mv "Solo.out/GeneFull/raw/matrix.mtx"   matrix.mtx
      else
        mv "Solo.out/GeneFull/raw/barcodes.tsv" barcodes.tsv
        mv "Solo.out/GeneFull/raw/features.tsv" features.tsv
        mv "Solo.out/GeneFull/raw/matrix.mtx"   matrix.mtx
        mv "Solo.out/Gene/raw/barcodes.tsv"     barcodes_sn_rna.tsv
        mv "Solo.out/Gene/raw/features.tsv"     features_sn_rna.tsv
        mv "Solo.out/Gene/raw/matrix.mtx"       matrix_sn_rna.mtx
      fi
    else
      echo Error: unknown counting mode: "$counting_mode". Should be either sn_rna or sc_rna.
    fi
    mv Aligned.sortedByCoord.out.bam ~{output_bam_basename}.bam

  >>>

  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = "~{output_bam_basename}.bam"
    File alignment_log = "Log.final.out"
    File general_log = "Log.out"
    File barcodes = "barcodes.tsv"
    File features = "features.tsv"
    File matrix = "matrix.mtx"
    File barcodes_sn_rna = "barcodes_sn_rna.tsv"
    File features_sn_rna = "features_sn_rna.tsv"
    File matrix_sn_rna = "matrix_sn_rna.mtx"
  }
}

task MergeStarOutput {

  input {
    Array[File] barcodes
    Array[File] features
    Array[File] matrix
    String input_id

    #runtime values
    String docker = "us.gcr.io/broad-gotc-prod/pytools:1.0.0-1661263730"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil(size(matrix, "Gi") * 2) + 10
    Int preemptible = 3
  }
  meta {
    description: "Create three files as .npy,  .npy and .npz for numpy array and scipy csr matrix for the barcodes, gene names and the count matrix by merging multiple  STARSolo count matrices (mtx format)."
  }

  parameter_meta {
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -e
    declare -a barcodes_files=(~{sep=' ' barcodes})
    declare -a features_files=(~{sep=' ' features})
    declare -a matrix_files=(~{sep=' ' matrix})

   # create the  compressed raw count matrix with the counts, gene names and the barcodes
    python3 /usr/gitc/create-merged-npz-output.py \
        --barcodes ${barcodes_files[@]} \
        --features ${features_files[@]} \
        --matrix ${matrix_files[@]} \
        --input_id ~{input_id}
  >>>

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File row_index = "~{input_id}_sparse_counts_row_index.npy"
    File col_index = "~{input_id}_sparse_counts_col_index.npy"
    File sparse_counts = "~{input_id}_sparse_counts.npz"
  }
}

task STARsoloFastqSlideSeq {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    File tar_star_reference
    File whitelist
    String output_bam_basename
    String read_structure
    Boolean? count_exons

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658781884"
    Int machine_mem_mb = 64000
    Int cpu = 8
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 3)) + ceil(size(r1_fastq, "Gi") * 20) +  ceil(size(r2_fastq, "Gi") * 20)
    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 3
  }

  command <<<
    set -e
    declare -a fastq1_files=(~{sep=' ' r1_fastq})
    declare -a fastq2_files=(~{sep=' ' r2_fastq})
    cut -f 1 ~{whitelist} > WhiteList.txt

    nums=$(echo ~{read_structure} | sed 's/[[:alpha:]]/ /g')
    read -a arr_num <<< $nums

    chars=$(echo ~{read_structure} | sed 's/[[:digit:]]/ /g')
    read -a arr_char <<< $chars

    UMILen=0
    CBLen=0
    for (( i=0; i<${#arr_char[@]}; ++i));
      do
        if [[ ${arr_char[$i]} == 'C' ]]
        then
          CBLen=$(( CBLen + arr_num[$i] ))
        elif [[ ${arr_char[$i]} == 'M' ]]
        then
          UMILen=$(( UMILen + arr_num[$i] ))
        fi
    done;
    UMIstart=$(( 1 + CBLen))

    # If this argument is true, we will count reads aligned to exons in addition
    COUNTING_MODE="GeneFull"
    if ~{count_exons}
    then
      COUNTING_MODE="Gene GeneFull"
    fi

    # prepare reference
    mkdir genome_reference
    tar -xf "~{tar_star_reference}" -C genome_reference --strip-components 1
    rm "~{tar_star_reference}"

    STAR \
      --soloType Droplet \
      --soloCBwhitelist WhiteList.txt \
      --soloFeatures $COUNTING_MODE \
      --runThreadN ~{cpu} \
      --genomeDir genome_reference \
      --readFilesIn $fastq2_files $fastq1_files \
      --readFilesCommand "gunzip -c" \
      --soloInputSAMattrBarcodeSeq CR UR \
      --soloInputSAMattrBarcodeQual CY UY \
      --soloCBlen $CBLen \
      --soloCBstart 1 \
      --soloUMIlen $UMILen \
      --soloUMIstart $UMIstart \
      --outSAMtype BAM SortedByCoordinate \
      --clip3pAdapterSeq AAAAAA \
      --clip3pAdapterMMp 0.1 \
      --outSAMattributes UB UR UY CR CB CY NH GX GN

    touch barcodes_exon.tsv
    touch features_exon.tsv
    touch matrix_exon.mtx

    mv "Solo.out/GeneFull/raw/barcodes.tsv" barcodes.tsv
    mv "Solo.out/GeneFull/raw/features.tsv" features.tsv
    mv "Solo.out/GeneFull/raw/matrix.mtx"   matrix.mtx

    if  ~{count_exons}
    then
      mv "Solo.out/Gene/raw/barcodes.tsv"     barcodes_exon.tsv
      mv "Solo.out/Gene/raw/features.tsv"     features_exon.tsv
      mv "Solo.out/Gene/raw/matrix.mtx"       matrix_exon.mtx
    fi

    mv Aligned.sortedByCoord.out.bam ~{output_bam_basename}.bam

  >>>

  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = "~{output_bam_basename}.bam"
    File alignment_log = "Log.final.out"
    File general_log = "Log.out"
    File barcodes = "barcodes.tsv"
    File features = "features.tsv"
    File matrix = "matrix.mtx"
    File barcodes_sn_rna = "barcodes_exon.tsv"
    File features_sn_rna = "features_exon.tsv"
    File matrix_sn_rna = "matrix_exon.mtx"  
  }
}

task STARGenomeRefVersion {
  input {
    String tar_star_reference
    Int disk = 10
  }

  meta {
    description: "Reads the reference file name and outputs a txt file containing genome source, build, and annotation version"
  }

  parameter_meta {
    tar_star_reference: "input STAR reference in TAR format"
  }

  command <<<
    # check genomic reference version and print to output txt file
    STRING=~{tar_star_reference}
    BASE=$(basename $STRING .tar)
    IFS=' -' read -r -a array <<< $BASE
    REFERENCE=${array[2]}
    VERSION=${array[4]}
    ANNOTATION=${array[5]}

    echo -e "$REFERENCE\n$VERSION\n$ANNOTATION" > reference_version.txt
    echo Reference is $REFERENCE
    echo Version is $VERSION
    echo Annotation is $ANNOTATION

  >>>

# Output is TXT file containing reference source, build version and annotation version
  output {
    File genomic_ref_version = "reference_version.txt"
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    memory: "2 GiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu:"1"
  }
}
