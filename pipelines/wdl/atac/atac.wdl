version 1.0

import "../../../tasks/wdl/MergeSortBam.wdl" as Merge
import "../../../tasks/wdl/FastqProcessing.wdl" as FastqProcessing
import "../../../tasks/wdl/PairedTagUtils.wdl" as AddBB
import "../../../tasks/wdl/Utilities.wdl" as utils
import "../peak_calling/PeakCalling.wdl" as peakcalling # import peakcalling as subworkflow

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
    # Additional library aliquot ID
    String? atac_nhash_id

    #Expected cells from library preparation. Also the number of barcodes kept in the filtered h5ad.
    Int atac_expected_cells = 3000
    # Anchor for SnapATAC2's OrdMag cell caller (number_of_cells metric only; does not filter)
    Int ordmag_anchor_cells = 10000
    # Write <input_id>.filtered.metrics.h5ad with the top atac_expected_cells barcodes
    Boolean filter_cells_to_expected = true

    # Option for running files with preindex
    Boolean preindex = false
    # Option for running peak calling, library level peak calling is always run 
    Boolean peak_calling = false
    
    # BWA ref
    File tar_bwa_reference
    # BWA machine type -- to select number of splits 
    Int num_threads_bwa = 128
    Int mem_size_bwa = 512
    String cpu_platform_bwa = "Intel Ice Lake"
    String vm_size

    # Text file containing chrom_sizes for genome build (i.e. hg38)
    File chrom_sizes
    #File for annotations for calculating ATAC TSSE
    File annotations_gtf
    # Whitelist
    File whitelist

    # TrimAdapters input
    String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"

    # Optional Aligned BAM input to skip alignment step
    File? aligned_ATAC_bam
  }

  String pipeline_version = "2.10.0"

  # Determine docker prefix based on cloud provider
  String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
  String acr_docker_prefix = "dsppipelinedev.azurecr.io/"
  String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

  # Docker image names
  String warp_tools_docker = "warp-tools:2.7.1"
  String cutadapt_docker = "cutadapt:1.0.0-4.4-1686752919"
  String samtools_docker = "samtools-dist-bwa:3.0.0"
  String upstools_docker = "upstools:1.0.0-2023.03.03-1704300311"
  String snap_atac_docker = "snapatac2:2.0.0"

  # Make sure either 'gcp' or 'azure' is supplied as cloud_provider input. If not, raise an error
  if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
    call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
        input:
            message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
    }
  }

  parameter_meta {
    read1_fastq_gzipped: "read 1 FASTQ file as input for the pipeline, contains read 1 of paired reads"
    read2_fastq_gzipped: "read 2 FASTQ file as input for the pipeline, contains the cellular barcodes corresponding to the reads in the read1 FASTQ and read 3 FASTQ"
    read3_fastq_gzipped: "read 3 FASTQ file as input for the pipeline, contains read 2 of paired reads"
    output_base_name: "base name to be used for the pipelines output and intermediate files"
    tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"
    num_threads_bwa: "Number of threads for bwa-mem2 task (default: 128)"
    mem_size_bwa: "Memory size in GB for bwa-mem2 task (default: 512)"
    cpu_platform_bwa: "CPU platform for bwa-mem2 task (default: Intel Ice Lake)"
    atac_expected_cells: "Number of cells loaded for library preparation. Used as the denominator of atac_percent_target and as the number of barcodes kept in the filtered h5ad (default: 3000)"
    ordmag_anchor_cells: "Anchor passed to SnapATAC2's OrdMag cell caller for the number_of_cells library metric; the reference barcode is the one at rank 1% of this value. Does not change which barcodes are written (default: 10000)"
    filter_cells_to_expected: "When true, CreateFragmentFile also writes a filtered h5ad with the top atac_expected_cells barcodes ranked by fragments overlapping peaks; the unfiltered per-barcode h5ad is always written (default: true)"
 }

  if (!defined(aligned_ATAC_bam)) {
    call GetNumSplits {
      input:
         nthreads = num_threads_bwa,
         mem_size = mem_size_bwa,
         cpu_platform = cpu_platform_bwa,
         vm_size = vm_size
    }

    call FastqProcessing.FastqProcessATAC as SplitFastq {
      input:
        read1_fastq = read1_fastq_gzipped,
        read3_fastq = read3_fastq_gzipped,
        barcodes_fastq = read2_fastq_gzipped,
        output_base_name = input_id,
        num_output_files = GetNumSplits.ranks_per_node_out,
        whitelist = whitelist,
        docker_path = docker_prefix + warp_tools_docker
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
          docker_path = docker_prefix + samtools_docker,
          cloud_provider = cloud_provider,
          vm_size = vm_size
    }
  }

  File aligned_bam = select_first([aligned_ATAC_bam, BWAPairedEndAlignment.bam_aligned_output])

  if (preindex) {
    call AddBB.AddBBTag as BBTag {
      input:
        bam = aligned_bam,
        input_id = input_id,
        docker_path = docker_prefix + upstools_docker
    }
    call CreateFragmentFile as BB_fragment {
      input:
        bam = BBTag.bb_bam,
        chrom_sizes = chrom_sizes,
        annotations_gtf = annotations_gtf,
        preindex = preindex,
        docker_path = docker_prefix + snap_atac_docker,
        atac_nhash_id = atac_nhash_id,
        atac_expected_cells = atac_expected_cells,
        ordmag_anchor_cells = ordmag_anchor_cells,
        filter_cells_to_expected = filter_cells_to_expected,
        input_id = input_id
    }
  }
  if (!preindex) {
    call CreateFragmentFile {
      input:
        bam = aligned_bam,
        chrom_sizes = chrom_sizes,
        annotations_gtf = annotations_gtf,
        preindex = preindex,
        docker_path = docker_prefix + snap_atac_docker,
        atac_nhash_id = atac_nhash_id,
        atac_expected_cells = atac_expected_cells,
        ordmag_anchor_cells = ordmag_anchor_cells,
        filter_cells_to_expected = filter_cells_to_expected,
        input_id = input_id
    }
    if (peak_calling) {
      # Peak calling runs on the filtered barcodes when filtering is enabled
      call peakcalling.PeakCalling as PeakCalling{
        input:
          output_base_name = input_id,
          annotations_gtf = annotations_gtf,
          metrics_h5ad = select_first([CreateFragmentFile.Snap_metrics_filtered, CreateFragmentFile.Snap_metrics]),
          chrom_sizes = chrom_sizes,
          cloud_provider = cloud_provider,
      }
    }
  }
  
  File bam_aligned_output_atac = select_first([BBTag.bb_bam, aligned_bam])
  File fragment_file_atac = select_first([BB_fragment.fragment_file, CreateFragmentFile.fragment_file])
  File fragment_file_index_atac = select_first([BB_fragment.fragment_file_index, CreateFragmentFile.fragment_file_index])
  File snap_metrics_atac = select_first([BB_fragment.Snap_metrics,CreateFragmentFile.Snap_metrics])
  File? snap_metrics_filtered_atac = if preindex then BB_fragment.Snap_metrics_filtered else CreateFragmentFile.Snap_metrics_filtered
  File library_metrics = select_first([BB_fragment.atac_library_metrics, CreateFragmentFile.atac_library_metrics])

  output {
    File bam_aligned_output = bam_aligned_output_atac
    File fragment_file = fragment_file_atac
    File fragment_file_index = fragment_file_index_atac
    # Unfiltered per-barcode h5ad (every barcode with a fragment; obs['is_cell'] and obs['is_retained'] flags)
    File snap_metrics = snap_metrics_atac
    # Top atac_expected_cells barcodes; only produced when filter_cells_to_expected is true
    File? snap_metrics_filtered = snap_metrics_filtered_atac
    File library_metrics_file = library_metrics
    File? cellbybin_h5ad_file = PeakCalling.cellbybin_h5ad
    File? cellbypeak_h5ad_file = PeakCalling.cellbypeak_h5ad
  }
}

# get number of splits
task GetNumSplits {
  input {
    # machine specs for bwa-mem2 task 
    Int nthreads
    Int mem_size
    String cpu_platform 
    String docker_image = "ubuntu@sha256:2e863c44b718727c860746568e1d54afd13b2fa71b160f5cd9058fc436217b30"
    String vm_size
  }

  parameter_meta {
    docker_image: "the ubuntu docker image (default: ubuntu@sha256:2e863c44b718727c860746568e1d54afd13b2fa71b160f5cd9058fc436217b30)"
    nthreads: "Number of threads per node (default: 128)"
    mem_size: "the size of memory used during alignment"
    vm_size: "the virtual machine used for the task"
  }

  command <<<
    set -euo pipefail
    echo "Get number of splits for bwa-mem2"
    echo "#############################################"
    echo "Machine specs for bwa-mem2 task"
    echo "#############################################"

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
    vm_size: vm_size
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
    docker_path: "The docker image path containing the runtime environment for this task"
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
    String reference_path = tar_bwa_reference
    String read_group_id = "RG1"
    String read_group_sample_name = "RGSN1"
    String suffix = "trimmed_adapters.fastq.gz"
    String output_base_name
    String docker_path
    String cloud_provider

    # Runtime attributes
    Int disk_size = 2000
    Int nthreads
    Int mem_size
    String cpu_platform
    String vm_size
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
    docker_path: "The docker image path containing the runtime environment for this task"
    cloud_provider: "The cloud provider for the pipeline."
    vm_size: "the virtual machine used for the task"
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
    echo "Reheading BAM with reference"
    /usr/temp/Open-Omics-Acceleration-Framework/applications/samtools/samtools view -H final.sorted.bam > header.txt
    echo -e "@CO\tReference genome used: ~{reference_path}" >> header.txt
    /usr/temp/Open-Omics-Acceleration-Framework/applications/samtools/samtools reheader header.txt final.sorted.bam > final.sorted.reheader.bam
    mv final.sorted.reheader.bam ~{bam_aligned_output_name}
        
    echo "the present working dir"
    pwd

    # save output logs for bwa-mem2
    mkdir output_logs
    mv *.txt output_logs

    if [ "~{cloud_provider}" == "gcp" ]; then
        tar -zcvf output_distbwa_log.tar.gz output_logs
        mv output_distbwa_log.tar.gz ../
    else
        tar -zcvf output_distbwa_log.tar.gz output_logs
        mv output_distbwa_log.tar.gz ../
    fi

    # move bam file to the root of cromwell
    # if the cloud provider is azure, move the file to /cromwell-executions
    # if the cloud provider is gcp, move the file to /cromwell_root
    if [ "~{cloud_provider}" == "gcp" ]; then
      mv ~{bam_aligned_output_name} ../
    else
      mv ~{bam_aligned_output_name} ../
    fi
  >>>

  runtime {
    docker: docker_path
    disks: "local-disk ${disk_size} SSD"
    cpu: nthreads
    cpuPlatform: cpu_platform
    memory: "${mem_size} GiB"
    vm_size: vm_size
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
    Array[String] mito_list = ['chrM', 'M']
    Int disk_size = 500
    Int mem_size = 64
    Int nthreads = 4
    String cpuPlatform = "Intel Cascade Lake"
    String docker_path
    String atac_nhash_id = ""
    String input_id
    Int atac_expected_cells = 3000
    Int ordmag_anchor_cells = 10000
    Boolean filter_cells_to_expected = true
    String gtf_path = annotations_gtf
  }

  parameter_meta {
    bam: "Aligned bam with CB in CB tag. This is the output of the BWAPairedEndAlignment task."
    chrom_sizes: "Text file containing chrom_sizes for genome build (i.e. hg38)."
    annotations_gtf: "GTF for SnapATAC2 to calculate TSS sites of fragment file."
    disk_size: "Disk size used in create fragment file step."
    mem_size: "The size of memory used in create fragment file."
    docker_path: "The docker image path containing the runtime environment for this task"
    atac_expected_cells: "Number of barcodes kept in the filtered h5ad (top barcodes ranked by fragments overlapping peaks) and the denominator of atac_percent_target. Default 3000."
    ordmag_anchor_cells: "Anchor (recovered_cells) passed to SnapATAC2's OrdMag cell caller when computing the number_of_cells metric; the reference barcode is the one at rank 1% of this value. Does not affect which barcodes are written. Default 10000."
    filter_cells_to_expected: "When true, also write <input_id>.filtered.metrics.h5ad containing the top atac_expected_cells barcodes. The unfiltered h5ad is always written. Default true."
  }

  command <<<
    set -euo pipefail
    set -x 

    python3 <<CODE

    # import libraries
    import snapatac2 as snap
    import numpy as np
    from collections import OrderedDict
    from snapatac2.preprocessing._cell_calling import filter_cellular_barcodes_ordmag
    import csv

    # set parameters
    bam = "~{bam}"
    input_id = "~{input_id}"
    atac_gtf = "~{annotations_gtf}"
    preindex = "~{preindex}"
    atac_nhash_id = "~{atac_nhash_id}"
    mito_list = "~{sep=' ' mito_list}"
    expected_cells = ~{atac_expected_cells}
    ordmag_anchor_cells = ~{ordmag_anchor_cells}
    filter_cells_to_expected = "~{filter_cells_to_expected}" == "true"
    gtf_path = "~{gtf_path}"

    mito_list = mito_list.split(" ")
    print("Mitochondrial chromosomes:", mito_list)

    # calculate chrom size dictionary based on text file
    chrom_size_dict = {}
    with open('~{chrom_sizes}', 'r') as f:
      for line in f:
        key, value = line.strip().split()
        chrom_size_dict[str(key)] = int(value)

    barcode_tag = "BB" if preindex == "true" else "CB"
    fragment_file = input_id + ".fragments.tsv"
    unfiltered_h5ad = input_id + ".metrics.h5ad"
    filtered_h5ad = input_id + ".filtered.metrics.h5ad"

    # The steps below reproduce snapatac2.pp.recipe_10x_metrics (SnapATAC2 2.7.1,
    # preprocessing/_recipe.py) with two deliberate differences:
    #   1. Cell calling is a METRIC only. The recipe's OrdMag call used
    #      recovered_cells=None, whose self-estimate collapses onto a few very deep
    #      barcodes; here it is anchored at ordmag_anchor_cells. The called set is
    #      flagged in obs['is_cell'] but is NOT used to subset the data.
    #   2. The recipe subsets the backed h5ad in place to the called cells, which
    #      truncated the per-barcode output. Here the unfiltered table is always
    #      written, and an optional filtered copy keeps the top atac_expected_cells
    #      barcodes ranked by fragments overlapping peaks.
    qc = OrderedDict()
    qc["NHashID"] = atac_nhash_id
    qc["Sequencing"] = {}
    qc["Cells"] = {}
    qc["Library Complexity"] = {}
    qc["Mapping"] = {}
    qc["Targeting"] = {}

    # 1. BAM -> fragment file (extract CB or BB tag; same defaults as the recipe)
    print("Creating fragment file from", barcode_tag, "tag")
    bam_qc = snap.pp.make_fragment_file(bam, fragment_file, is_paired=True, barcode_tag=barcode_tag, chrM=mito_list)
    qc["Sequencing"]["Sequenced_reads"] = bam_qc["sequenced_reads"]
    qc["Sequencing"]["Sequenced_read_pairs"] = bam_qc["sequenced_read_pairs"]
    qc["Sequencing"]["Fraction_valid_barcode"] = bam_qc["frac_valid_barcode"]
    qc["Sequencing"]["Fraction_Q30_bases_in_read_1"] = bam_qc["frac_q30_bases_read1"]
    qc["Sequencing"]["Fraction_Q30_bases_in_read_2"] = bam_qc["frac_q30_bases_read2"]
    qc["Mapping"]["Fraction_confidently_mapped"] = bam_qc["frac_confidently_mapped"]
    qc["Mapping"]["Fraction_unmapped"] = bam_qc["frac_unmapped"]
    qc["Mapping"]["Fraction_nonnuclear"] = bam_qc["frac_nonnuclear"]
    qc["Mapping"]["Fraction_fragment_in_nucleosome_free_region"] = bam_qc["frac_fragment_in_nucleosome_free_region"]
    qc["Mapping"]["Fraction_fragment_flanking_single_nucleosome"] = bam_qc["frac_fragment_flanking_single_nucleosome"]
    qc["Library Complexity"]["Fraction_duplicates"] = bam_qc["frac_duplicates"]

    # 2. Fragment file -> per-barcode table, backed by the unfiltered h5ad.
    #    min_num_fragments=0 keeps every barcode with at least one fragment.
    print("Importing fragments into", unfiltered_h5ad)
    adata = snap.pp.import_data(fragment_file, chrom_sizes=chrom_size_dict, file=unfiltered_h5ad, min_num_fragments=0, chrM=mito_list)
    n_barcodes = adata.n_obs
    print("Barcodes imported:", n_barcodes)

    # 3. TSS enrichment (per barcode in obs['tsse'] and library level in uns), excluding mito contigs
    snap.metrics.tsse(adata, atac_gtf, exclude_chroms=mito_list)
    qc["Targeting"]["TSS_enrichment_score"] = adata.uns['library_tsse']
    qc["Targeting"]["Fraction_of_high-quality_fragments_overlapping_TSS"] = adata.uns['frac_overlap_TSS']

    # 4. Peaks: MACS3 on the pseudobulk of all barcodes (recipe default qvalue)
    snap.tl.macs3(adata, qvalue=0.001)
    peaks = [f"{row[0]}:{row[1]}-{row[2]}" for row in adata.uns['macs3_pseudobulk'].iter_rows()]
    qc["Targeting"]["Number_of_peaks"] = len(peaks)
    qc["Targeting"]["Fraction_of_genome_in_peaks"] = snap._snapatac2.total_size_of_peaks(peaks) / adata.uns['reference_sequences']['reference_seq_length'].sum()

    # 5. Fragments overlapping peaks per barcode (stored in obs['n_frag_overlap_peak'])
    snap.metrics.frip(adata, {"n_frag_overlap_peak": peaks}, normalized=False)
    n_fragment = adata.obs['n_fragment'].to_numpy()
    n_frag_in_peaks = adata.obs['n_frag_overlap_peak'].to_numpy()
    qc["Targeting"]["Fraction_of_high-quality_fragments_overlapping_peaks"] = n_frag_in_peaks.sum() / n_fragment.sum()

    # 6. Cell calling for the library metric: OrdMag anchored at ordmag_anchor_cells.
    #    This is snapatac2.pp.call_cells with recovered_cells set instead of None.
    cell_idx, ordmag_stats = filter_cellular_barcodes_ordmag(n_frag_in_peaks.astype(int), ordmag_anchor_cells)
    cell_idx = np.sort(np.asarray(cell_idx, dtype=int))
    n_cells = int(len(cell_idx))
    print("OrdMag called", n_cells, "cells with anchor", ordmag_anchor_cells)
    is_paired = bam_qc["sequenced_read_pairs"] > 0
    qc["Cells"]["Number_of_cells"] = n_cells
    qc["Cells"]["Ordmag_anchor_cells"] = ordmag_anchor_cells
    if n_cells > 0:
        raw = bam_qc["sequenced_read_pairs"] if is_paired else bam_qc["sequenced_reads"]
        qc["Cells"]["Mean_raw_read_pairs_per_cell"] = raw / n_cells
        qc["Cells"]["Median_high-quality_fragments_per_cell"] = float(np.median(n_fragment[cell_idx]))
        qc["Cells"]["Fraction of high-quality fragments in cells"] = n_fragment[cell_idx].sum() / n_fragment.sum()
        qc["Cells"]["Min_fragments_in_peaks_in_cells"] = int(n_frag_in_peaks[cell_idx].min())
        # transposition events (both fragment ends) in peaks, summed over called cells
        insertions = snap.metrics.frip(adata, {"overlap_peak": peaks}, normalized=False, count_as_insertion=True, inplace=False)
        insertions_in_peaks = np.asarray(insertions["overlap_peak"])
        denominator = n_fragment[cell_idx].sum() * (2 if is_paired else 1)
        qc["Cells"]["Fraction_of_transposition_events_in_peaks_in_cells"] = insertions_in_peaks[cell_idx].sum() / denominator
    else:
        for key in ["Mean_raw_read_pairs_per_cell", "Median_high-quality_fragments_per_cell",
                    "Fraction of high-quality fragments in cells", "Min_fragments_in_peaks_in_cells",
                    "Fraction_of_transposition_events_in_peaks_in_cells"]:
            qc["Cells"][key] = float("nan")
    qc["Cells"]["atac_percent_target"] = n_cells / expected_cells * 100

    # 7. Annotate the per-barcode table. is_cell = OrdMag call; is_retained = filter below.
    is_cell = np.zeros(n_barcodes, dtype=bool)
    is_cell[cell_idx] = True
    adata.obs['is_cell'] = is_cell
    adata.uns['NHashID'] = atac_nhash_id
    adata.uns['reference_gtf_file'] = gtf_path

    # 8. Filtering, decoupled from cell calling: keep the top atac_expected_cells
    #    barcodes by fragments overlapping peaks (ties broken by import order).
    n_retained = min(int(expected_cells), n_barcodes) if filter_cells_to_expected else 0
    is_retained = np.zeros(n_barcodes, dtype=bool)
    if n_retained > 0:
        order = np.argsort(-n_frag_in_peaks, kind="stable")
        keep_idx = np.sort(order[:n_retained])
        is_retained[keep_idx] = True
        qc["Cells"]["Number_of_barcodes_retained"] = int(n_retained)
        qc["Cells"]["Min_fragments_in_peaks_retained"] = int(n_frag_in_peaks[keep_idx].min())
    adata.obs['is_retained'] = is_retained

    if n_retained > 0:
        print("Writing filtered h5ad with", n_retained, "barcodes to", filtered_h5ad)
        filtered = adata.subset(keep_idx, out=filtered_h5ad, inplace=False)
        filtered.close()
    adata.close()

    # Flatten the metrics dictionary into (metric, value) rows
    flattened_data = []
    for category, metrics in qc.items():
        if isinstance(metrics, dict):
            for metric, value in metrics.items():
                flattened_data.append((metric, value))
        else:
            flattened_data.append((category, metrics))

    # Convert the flattened keys to lowercase (except for 'NHashID')
    flattened_data = [(metric if metric == 'NHashID' else str(metric).lower(), value) for metric, value in flattened_data]

    # Write to CSV
    csv_file_path = "~{input_id}_~{atac_nhash_id}_library_metrics.csv"
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(flattened_data)

    print(f"Dictionary successfully written to {csv_file_path}")

    CODE
    
    # sorting the file
    echo "Sorting file"
    sort -k1,1V -k2,2n "~{input_id}.fragments.tsv" > "~{input_id}.fragments.sorted.tsv"
    echo "Starting bgzip"
    bgzip "~{input_id}.fragments.sorted.tsv"
    echo "Starting tabix"
    tabix -s 1 -b 2 -e 3 -C "~{input_id}.fragments.sorted.tsv.gz"
  >>>

  runtime {
    docker: docker_path
    disks: "local-disk ${disk_size} SSD"
    memory: "${mem_size} GiB"
    cpu: nthreads
    cpuPlatform: cpuPlatform
  }

  output {
    File fragment_file = "~{input_id}.fragments.sorted.tsv.gz"
    File fragment_file_index = "~{input_id}.fragments.sorted.tsv.gz.csi"
    File Snap_metrics = "~{input_id}.metrics.h5ad"
    File? Snap_metrics_filtered = "~{input_id}.filtered.metrics.h5ad"
    File atac_library_metrics = "~{input_id}_~{atac_nhash_id}_library_metrics.csv"
  }
}
