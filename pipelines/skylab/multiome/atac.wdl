version 1.0

import "../../../tasks/skylab/MergeSortBam.wdl" as Merge
import "../../../tasks/skylab/FastqProcessing.wdl" as FastqProcessing
import "../../../tasks/skylab/PairedTagUtils.wdl" as AddBB

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
    String cloud_provider

    # Option for running files with preindex
    Boolean preindex = false
    
    # BWA ref
    File tar_bwa_reference
    # BWA machine type -- to select number of splits 
    Int num_threads_bwa = 128
    Int mem_size_bwa = 512
    String cpu_platform_bwa = "Intel Ice Lake"

    # GTF for SnapATAC2 to calculate TSS sites of fragment file
    File annotations_gtf
    # Text file containing chrom_sizes for genome build (i.e. hg38)
    File chrom_sizes
    # Whitelist
    File whitelist

    # TrimAdapters input
    String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
  }

  String pipeline_version = "1.1.7"

  # Determine docker prefix based on cloud provider
  String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
  String acr_docker_prefix = "dsppipelinedev.azurecr.io/"
  String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

  # Docker image names
  String warp_tools_2_0_0 = "warp-tools:2.0.0"
  String cutadapt_docker = "cutadapt:1.0.0-4.4-1686752919"
  String sam_tools_docker = "samtools-dist-bwa:2.0.0"
  String upstools_docker = "upstools:1.0.0-2023.03.03-1704300311"
  String snap_atac_docker = "snapatac2:1.0.4-2.3.1"

  parameter_meta {
    read1_fastq_gzipped: "read 1 FASTQ file as input for the pipeline, contains read 1 of paired reads"
    read2_fastq_gzipped: "read 2 FASTQ file as input for the pipeline, contains the cellular barcodes corresponding to the reads in the read1 FASTQ and read 3 FASTQ"
    read3_fastq_gzipped: "read 3 FASTQ file as input for the pipeline, contains read 2 of paired reads"
    output_base_name: "base name to be used for the pipelines output and intermediate files"
    tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"
    num_threads_bwa: "Number of threads for bwa-mem2 task (default: 128)"
    mem_size_bwa: "Memory size in GB for bwa-mem2 task (default: 512)"
    cpu_platform_bwa: "CPU platform for bwa-mem2 task (default: Intel Ice Lake)"
  
 }

  call GetNumSplits {
    input:
       nthreads = num_threads_bwa, 
       mem_size = mem_size_bwa,
       cpu_platform = cpu_platform_bwa
  }

  call FastqProcessing.FastqProcessATAC as SplitFastq {
    input:
      read1_fastq = read1_fastq_gzipped,
      read3_fastq = read3_fastq_gzipped,
      barcodes_fastq = read2_fastq_gzipped,
      output_base_name = input_id,
      num_output_files = GetNumSplits.ranks_per_node_out,
      whitelist = whitelist,
      docker_path = docker_prefix + warp_tools_2_0_0
  }

  scatter(idx in range(length(SplitFastq.fastq_R1_output_array))) {
    call TrimAdapters {
      input:
        read1_fastq = SplitFastq.fastq_R1_output_array[idx],
        read3_fastq = SplitFastq.fastq_R3_output_array[idx],
        output_base_name = input_id + "_" + idx,
        adapter_seq_read1 = adapter_seq_read1,
        adapter_seq_read3 = adapter_seq_read3,
        docker_path = docker_prefix + cutadapt_docker
    }
  }

  call BWAPairedEndAlignment {
    input:
        read1_fastq = TrimAdapters.fastq_trimmed_adapter_output_read1,
        read3_fastq = TrimAdapters.fastq_trimmed_adapter_output_read3,
        tar_bwa_reference = tar_bwa_reference,
        output_base_name = input_id,
        nthreads = num_threads_bwa, 
        mem_size = mem_size_bwa,
        cpu_platform = cpu_platform_bwa,
        docker_path = docker_prefix + sam_tools_docker
  }

  if (preindex) {
    call AddBB.AddBBTag as BBTag {
      input:
        bam = BWAPairedEndAlignment.bam_aligned_output,
        input_id = input_id,
        docker_path = docker_prefix + upstools_docker
    }
    call CreateFragmentFile as BB_fragment {
      input:
        bam = BBTag.bb_bam,
        chrom_sizes = chrom_sizes,
        annotations_gtf = annotations_gtf,
        preindex = preindex,
        docker_path = docker_prefix + snap_atac_docker
    }
  }
  if (!preindex) {
    call CreateFragmentFile {
      input:
        bam = BWAPairedEndAlignment.bam_aligned_output,
        chrom_sizes = chrom_sizes,
        annotations_gtf = annotations_gtf,
        preindex = preindex,
        docker_path = docker_prefix + snap_atac_docker

    }
  }
  File bam_aligned_output_atac = select_first([BBTag.bb_bam, BWAPairedEndAlignment.bam_aligned_output])
  File fragment_file_atac = select_first([BB_fragment.fragment_file, CreateFragmentFile.fragment_file])
  File snap_metrics_atac = select_first([BB_fragment.Snap_metrics,CreateFragmentFile.Snap_metrics])

  output {
    File bam_aligned_output = bam_aligned_output_atac
    File fragment_file = fragment_file_atac
    File snap_metrics = snap_metrics_atac
  }
}

# get number of splits
task GetNumSplits {
  input {
    # machine specs for bwa-mem2 task 
    Int nthreads
    Int mem_size
    String cpu_platform 
    String docker_image = "ubuntu:latest"
  }

  parameter_meta {
    docker_image: "the ubuntu docker image (default: ubuntu:latest)"
    nthreads: "Number of threads per node (default: 128)"
    mem_size: "the size of memory used during alignment"
  }

  command <<<
    set -euo pipefail
    
    # steps taken from https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/blob/main/pipelines/fq2sortedbam/print_config.sh
    num_nodes=1
    lscpu
    lscpu > compute_config
    
    num_cpus_per_node=$(cat compute_config | grep -E '^CPU\(s\)' | awk  '{print $2}')
    num_sockets=$(cat compute_config | grep -E '^Socket'| awk  '{print $2}')
    num_numa=$(cat compute_config | grep '^NUMA node(s)' | awk '{print $3}')
    num_cpus_all_node=`expr ${num_cpus_per_node} \* ${num_nodes}`
    threads_per_core=$(cat compute_config | grep -E '^Thread' | awk  '{print $4}')
    
    num_cpus_all_node=`expr ${num_cpus_per_node} \* ${num_nodes}`

    echo "Number of threads: " $num_cpus_per_node
    echo "Number of sockets: " $num_sockets
    echo "Number of NUMA domains: "$num_numa
    echo "Number of threads per core: "$threads_per_core
    echo "Number of CPUs: $num_cpus_all_node"

    num_physical_cores_all_nodes=`expr ${num_cpus_all_node} / ${threads_per_core}`
    num_physical_cores_per_nodes=`expr  ${num_cpus_per_node} / ${threads_per_core}`
    num_physical_cores_per_socket=`expr ${num_physical_cores_all_nodes} / ${num_sockets}`
    num_physical_cores_per_numa=`expr ${num_physical_cores_all_nodes} / ${num_numa}`
    echo "Number physical cores: "$num_physical_cores_per_nodes
    echo "Number physical cores per socket: "$num_physical_cores_per_socket
    echo "Number physical cores per numa: "$num_physical_cores_per_numa

    th=`expr ${num_physical_cores_per_numa} / 2` 
    if [ $th -le 10 ]
    then
        th=${num_physical_cores_per_numa}
    fi

    while [ $num_physical_cores_per_nodes -gt $th ]
    do
        num_physical_cores_per_nodes=`expr $num_physical_cores_per_nodes / 2`
    done

    num_physical_cores_per_rank=$num_physical_cores_per_nodes
    total_num_ranks=`expr ${num_physical_cores_all_nodes} / ${num_physical_cores_per_rank}`

    ranks_per_node=`expr ${total_num_ranks} / ${num_nodes}`
    echo "Number of MPI ranks: "${total_num_ranks}
    echo "Number of cores per MPI rank: "$num_physical_cores_per_nodes
    echo "#############################################"
    #echo "Note: Each MPI rank runs a bwa-mem2 process on its input fastq files produced by fqprocess. Please ensure that the number of files created due to bam_size parameter to fqprocess (in config file) creates number of fastq files equal to ${total_num_ranks}"
    echo "Please set bam_size such that fastqprocess creates ${total_num_ranks} splits of input fastq files"
    echo "#############################################"

    echo $total_num_ranks > total_num_ranks.txt
    echo $ranks_per_node > ranks_per_node.txt

  >>>

  runtime {
    docker: docker_image
    cpu: nthreads
    cpuPlatform: cpu_platform
    memory: "${mem_size} GiB"
  }

  output {
    Int ranks_per_node_out = read_int("ranks_per_node.txt")
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
    String docker_path
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
    docker: docker_path
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
    String docker_path

    # Runtime attributes
    Int disk_size = 2000
    Int nthreads
    Int mem_size
    String cpu_platform 
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

  String bam_aligned_output_name = output_base_name + ".bam"

  # bwa and call samtools to convert sam to bam
  command <<<

    set -euo pipefail    

    # print lscpu  
    echo "lscpu output"
    lscpu
    echo "end of lscpu output"

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
    
    PARAMS="+R '@RG\tID:~{read_group_id}\tSM:~{read_group_sample_name}' +C"
    
    INPUT_DIR=$file_path/input_dir
    OUTPUT_DIR=$file_path/output_dir
    
    input_to_config="INPUT_DIR=\"${INPUT_DIR}\"\nOUTPUT_DIR=\"${OUTPUT_DIR}\"\nPREFIX=\"${PREFIX}\"\nSUFFIX=\"${SUFFIX}\"\n"
    other_to_add="R1=\"${R1}\"\nR2=\"${R2}\"\nR3=\"${R3}\"\nI1=\"${I1}\"\nREF=\"${REF}\"\n"
    params="PARAMS=\"${PARAMS}\""
    
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
  >>>

  runtime {
    docker: docker_path
    disks: "local-disk ${disk_size} SSD"
    cpu: nthreads
    cpuPlatform: cpu_platform
    memory: "${mem_size} GiB"
  }

  output {
    File bam_aligned_output = bam_aligned_output_name
    File output_distbwa_log_tar = "output_distbwa_log.tar.gz"
  }
}

# make fragment file
task CreateFragmentFile {
  input {
    File bam
    File annotations_gtf
    File chrom_sizes
    Boolean preindex
    Int disk_size = 500
    Int mem_size = 16
    Int nthreads = 1
    String cpuPlatform = "Intel Cascade Lake"
    String docker_path
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
    preindex = "~{preindex}"

    # calculate chrom size dictionary based on text file
    chrom_size_dict={}
    with open('~{chrom_sizes}', 'r') as f:
      for line in f:
        key, value = line.strip().split()
        chrom_size_dict[str(key)] = int(value)

    # use snap atac2
    import snapatac2.preprocessing as pp
    import snapatac2 as snap

    # extract CB or BB (if preindex is true) tag from bam file to create fragment file
    if preindex == "true":
      pp.make_fragment_file("~{bam}", "~{bam_base_name}.fragments.tsv", is_paired=True, barcode_tag="BB")
    elif preindex == "false":
      pp.make_fragment_file("~{bam}", "~{bam_base_name}.fragments.tsv", is_paired=True, barcode_tag="CB")
      

    # calculate quality metrics; note min_num_fragments and min_tsse are set to 0 instead of default
    # those settings allow us to retain all barcodes
    pp.import_data("~{bam_base_name}.fragments.tsv", file="~{bam_base_name}.metrics.h5ad", chrom_size=chrom_size_dict, gene_anno="~{annotations_gtf}", min_num_fragments=0, min_tsse=0)

    CODE
  >>>

  runtime {
    docker: docker_path
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
