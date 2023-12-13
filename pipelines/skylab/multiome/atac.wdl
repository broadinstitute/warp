version 1.0

import "../../../tasks/skylab/MergeSortBam.wdl" as Merge
import "../../../tasks/skylab/FastqProcessing.wdl" as FastqProcessing

workflow ATAC {
  meta {
    description: "Processing for single-cell ATAC-seq data from the level of raw fastq reads. This is the first step of the multiome pipeline. ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a technique used in molecular biology to assess genome-wide chromatin accessibility. This pipeline processes 10x Genomics Multiome ATAC FASTQ files."
    allowNestedInputs: true
  }

  input {
    # Fastq inputs
    Array[String] read1_fastq_gzipped
    Array[String] read2_fastq_gzipped
    Array[String] read3_fastq_gzipped

    # Output prefix/base name for all intermediate files and pipeline outputs
    String input_id

    # BWA ref
    File tar_bwa_reference

    # GTF for SnapATAC2 to calculate TSS sites of fragment file
    File annotations_gtf
    # Text file containing chrom_sizes for genome build (i.e. hg38)
    File chrom_sizes
    # Number of fastq files for fastqprocess
    Int num_output_files = 4
    # Whitelist
    File whitelist

    # TrimAdapters input
    String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
  }

  String pipeline_version = "1.1.2"

  parameter_meta {
    read1_fastq_gzipped: "read 1 FASTQ file as input for the pipeline, contains read 1 of paired reads"
    read2_fastq_gzipped: "read 2 FASTQ file as input for the pipeline, contains the cellular barcodes corresponding to the reads in the read1 FASTQ and read 3 FASTQ"
    read3_fastq_gzipped: "read 3 FASTQ file as input for the pipeline, contains read 2 of paired reads"
    output_base_name: "base name to be used for the pipelines output and intermediate files"
    tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"
  }

  call FastqProcessing.FastqProcessATAC as SplitFastq {
    input:
      read1_fastq = read1_fastq_gzipped,
      read3_fastq = read3_fastq_gzipped,
      barcodes_fastq = read2_fastq_gzipped,
      output_base_name = input_id,
      num_output_files = num_output_files,
      whitelist = whitelist
  }

  scatter(idx in range(length(SplitFastq.fastq_R1_output_array))) {

    call TrimAdapters {
      input:
        read1_fastq = SplitFastq.fastq_R1_output_array[idx],
        read3_fastq = SplitFastq.fastq_R3_output_array[idx],
        output_base_name = input_id + "_" + idx,
        adapter_seq_read1 = adapter_seq_read1,
        adapter_seq_read3 = adapter_seq_read3
    }
  }

  call BWAPairedEndAlignment {
    input:
        read1_fastq = TrimAdapters.fastq_trimmed_adapter_output_read1,
        read3_fastq = TrimAdapters.fastq_trimmed_adapter_output_read3,
        tar_bwa_reference = tar_bwa_reference,
        output_base_name = input_id
  }

  call CreateFragmentFile {
    input:
      bam = BWAPairedEndAlignment.bam_aligned_output,
      chrom_sizes = chrom_sizes,
      annotations_gtf = annotations_gtf
  }

  output {
    File bam_aligned_output = BWAPairedEndAlignment.bam_aligned_output
    File fragment_file = CreateFragmentFile.fragment_file
    File snap_metrics = CreateFragmentFile.Snap_metrics
  }
}

# trim read 1 and read 2 adapter sequeunce with cutadapt
task TrimAdapters {
  input {
    File read1_fastq
    File read3_fastq
    String output_base_name

    Int min_length = 10
    Int quality_cutoff = 0

    String adapter_seq_read1
    String adapter_seq_read3

    # Runtime attributes/docker
    Int disk_size = ceil(2 * ( size(read1_fastq, "GiB") + size(read3_fastq, "GiB") )) + 200
    Int mem_size = 4
    String docker_image = "us.gcr.io/broad-gotc-prod/cutadapt:1.0.0-4.4-1686752919"
  }

  parameter_meta {
    read1_fastq: "read 1 fastq file containing sequencing reads as input for the pipeline"
    read3_fastq: "read 3 fastq file containing sequencing reads as input for the pipeline"
    min_length: "the minimum length for trimming. Reads that are too short even before adapter removal are also discarded"
    quality_cutoff: "cutadapt option to trim low-quality ends from reads before adapter removal"
    adapter_seq_read1: "cutadapt option for the sequence adapter for read 1 fastq"
    adapter_seq_read3: "cutadapt option for the sequence adapter for read 3 fastq"
    output_base_name: "base name to be used for the output of the task"
    docker_image: "the docker image using cutadapt to be used (default:us.gcr.io/broad-gotc-prod/cutadapt:1.0.0-4.4-1686752919)"
    mem_size: "the size of memory used during trimming adapters"
    disk_size : "disk size used in trimming adapters step"
  }

  # output names for trimmed reads
  String fastq_trimmed_adapter_output_name_read1 = output_base_name + ".R1.trimmed_adapters.fastq.gz"
  String fastq_trimmed_adapter_output_name_read3 = output_base_name + ".R3.trimmed_adapters.fastq.gz"

  # using cutadapt to trim off sequence adapters
  command <<<
    set -euo pipefail

    # fastq's, "-f", -A for paired adapters read 2"
    cutadapt \
    -Z \
    --minimum-length ~{min_length} \
    --quality-cutoff ~{quality_cutoff} \
    --adapter ~{adapter_seq_read1} \
    -A ~{adapter_seq_read3} \
    --output ~{fastq_trimmed_adapter_output_name_read1} \
    --paired-output ~{fastq_trimmed_adapter_output_name_read3} \
    ~{read1_fastq} ~{read3_fastq}
  >>>

  # use docker image for given tool cutadapat
  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    memory: "${mem_size} GiB"
  }

  output {
    File fastq_trimmed_adapter_output_read1 = fastq_trimmed_adapter_output_name_read1
    File fastq_trimmed_adapter_output_read3 = fastq_trimmed_adapter_output_name_read3
  }
}

# align the two trimmed fastq as paired end data using BWA
task BWAPairedEndAlignment {
  input {
    Array[File] read1_fastq
    Array[File] read3_fastq
    File tar_bwa_reference
    String read_group_id = "RG1"
    String read_group_sample_name = "RGSN1"
    String suffix = "trimmed_adapters.fastq.gz"
    String output_base_name
    String docker_image = "us.gcr.io/broad-gotc-prod/samtools-dist-bwa:aa-dist-bwa-"
    File monitoring = "gs://fc-51792410-8543-49ba-ad3f-9e274900879f/cromwell_monitoring_script2.sh"

    # Runtime attributes
    Int disk_size = 2000
    Int nthreads = 128
    Int mem_size = 512
  }

  parameter_meta {
    read1_fastq: "the trimmed read 1 fastq file containing sequencing reads as input for the aligner"
    read3_fastq: "the trimmed read 1 fastq file containing sequencing reads as input for the aligner"
    tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"
    read_group_id: "the read group id to be added upon alignment"
    read_group_sample_name: "the read group sample to be added upon alignment"
    nthreads: "the number of threads to use during bwa alignment"
    mem_size: "the size of memory used during alignment"
    disk_size : "disk size used in bwa alignment step"
    output_base_name: "basename to be used for the output of the task"
    docker_image: "the docker image using BWA to be used (default: us.gcr.io/broad-gotc-prod/samtools-bwa-mem-2:1.0.0-2.2.1_x64-linux-1685469504)"
  }

  String bam_aligned_output_name = output_base_name + ".aligned.bam"

  # bwa and call samtools to convert sam to bam
  command <<<

    set -euo pipefail    

    # print lscpu  
    echo "lscpu output"
    lscpu
    echo "end of lscpu output"

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring}" ]; then
      chmod a+x ~{monitoring}
      ~{monitoring} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # prepare reference
    declare -r REF_DIR=$(mktemp -d genome_referenceXXXXXX)
    tar -xf "~{tar_bwa_reference}" -C $REF_DIR --strip-components 1
    rm "~{tar_bwa_reference}"
    REF_PAR_DIR=$(basename "$(dirname "$REF_DIR/genome.fa")")
    echo $REF_PAR_DIR

    # make read1_fastq and read3_fastq into arrays 
    declare -a R1_ARRAY=(~{sep=' ' read1_fastq})
    declare -a R3_ARRAY=(~{sep=' ' read3_fastq})
    
    file_path=`pwd`
    echo "The current working directory is" $file_path

    # make input and output directories needed for distributed bwamem2 code
    mkdir "output_dir"
    mkdir "input_dir"
    
    echo "Move R1, R3 and reference files to input directory."
    R1=""
    echo "R1"
    for fastq in "${R1_ARRAY[@]}"; do mv "$fastq" input_dir; R1+=`basename $fastq`" "; done
    echo $R1
    R3=""
    echo "R3"
    for fastq in "${R3_ARRAY[@]}"; do mv "$fastq" input_dir; R3+=`basename $fastq`" "; done
    echo $R3

    mv $REF_DIR input_dir

    echo "List of files in input directory"
    ls input_dir
    
    # multiome-practice-may15_arcgtf, trimmed_adapters.fastq.gz
    PREFIX=~{output_base_name}
    SUFFIX=~{suffix}

    I1=""
    R2=""
    
    echo "REF_PAR_DIR:" $REF_PAR_DIR
    REF=$REF_PAR_DIR/genome.fa
    
    PARAMS='+R "@RG\tID:~{read_group_id}\tSM:~{read_group_sample_name}" +C'
    
    INPUT_DIR=$file_path/input_dir
    OUTPUT_DIR=$file_path/output_dir
    
    input_to_config="INPUT_DIR=\"${INPUT_DIR}\"\nOUTPUT_DIR=\"${OUTPUT_DIR}\"\nPREFIX=\"${PREFIX}\"\nSUFFIX=\"${SUFFIX}\"\n"
    other_to_add="R1=\"${R1}\"\nR2=\"${R2}\"\nR3=\"${R3}\"\nI1=\"${I1}\"\nREF=\"${REF}\"\n"
    params="PARAMS=\'${PARAMS}\'"
    
    printf "%b" "$input_to_config"
    printf "%b" "$other_to_add"
    echo $params
    
    # cd into fq2sortedbam
    cd /usr/temp/Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam
    # remove the first part of config
    tail -10 config > config
    # add inputs to config file (this file is needed to run bwa-mem2 in this specific code"
    printf "%b" "$input_to_config" | tee -a config
    printf "%b" "$other_to_add" | tee -a config
    echo $params | tee -a config
    echo "CONFIG"
    cat config
    # run bwa-mem2
    echo "Run distributed BWA-MEM2"
    ./run_bwa.sh multifq
    echo "Done running distributed BWA-MEM2"
    echo "List of files in output directory"
    ls $OUTPUT_DIR
    cd $OUTPUT_DIR
    
    # remove all files except for final and text file 
    echo "Remove all files except for final bam file and log files"
    ls | grep -xv final.sorted.bam | grep -v .txt$ | xargs rm

    echo "List of files in output directory after removal"
    ls
    
    # rename file to this
    mv final.sorted.bam ~{bam_aligned_output_name}
        
    # save output logs for bwa-mem2
    mkdir output_logs
    mv *txt output_logs
    tar -zcvf /cromwell_root/output_distbwa_log.tar.gz output_logs  
    
    # move bam file to /cromwell_root
    mv ~{bam_aligned_output_name} /cromwell_root
    ls /cromwell_root

  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} SSD"
    cpu: nthreads
    cpuPlatform: "Intel Ice Lake"
    memory: "${mem_size} GiB"
  }

  output {
    File bam_aligned_output = bam_aligned_output_name
    File monitoring_out = "monitoring.log"
    File output_distbwa_log_tar = "output_distbwa_log.tar.gz"
  }
}

# make fragment file
task CreateFragmentFile {
  input {
    File bam
    File annotations_gtf
    File chrom_sizes
    Int disk_size = 500
    Int mem_size = 16
    Int nthreads = 1
    String cpuPlatform = "Intel Cascade Lake"
  }

  String bam_base_name = basename(bam, ".bam")

  parameter_meta {
    bam: "Aligned bam with CB in CB tag. This is the output of the BWAPairedEndAlignment task."
    annotations_gtf: "GTF for SnapATAC2 to calculate TSS sites of fragment file."
    chrom_sizes: "Text file containing chrom_sizes for genome build (i.e. hg38)."
    disk_size: "Disk size used in create fragment file step."
    mem_size: "The size of memory used in create fragment file."
  }

  command <<<
    set -e pipefail

    python3 <<CODE

    # set parameters
    atac_gtf = "~{annotations_gtf}"
    bam = "~{bam}"
    bam_base_name = "~{bam_base_name}"
    chrom_sizes = "~{chrom_sizes}"

    # calculate chrom size dictionary based on text file
    chrom_size_dict={}
    with open('~{chrom_sizes}', 'r') as f:
      for line in f:
        key, value = line.strip().split()
        chrom_size_dict[str(key)] = int(value)

    # use snap atac2
    import snapatac2.preprocessing as pp
    import snapatac2 as snap

    # extract CB tag from bam file to create fragment file
    pp.make_fragment_file("~{bam}", "~{bam_base_name}.fragments.tsv", is_paired=True, barcode_tag="CB")

    # calculate quality metrics; note min_num_fragments and min_tsse are set to 0 instead of default
    # those settings allow us to retain all barcodes
    pp.import_data("~{bam_base_name}.fragments.tsv", file="~{bam_base_name}.metrics.h5ad", chrom_size=chrom_size_dict, gene_anno="~{annotations_gtf}", min_num_fragments=0, min_tsse=0)

    CODE
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/snapatac2:1.0.4-2.3.1"
    disks: "local-disk ${disk_size} SSD"
    memory: "${mem_size} GiB"
    cpu: nthreads
    cpuPlatform: cpuPlatform
  }

  output {
    File fragment_file = "~{bam_base_name}.fragments.tsv"
    File Snap_metrics = "~{bam_base_name}.metrics.h5ad"
  }
}
