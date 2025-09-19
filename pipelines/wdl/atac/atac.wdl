version 1.0

import "../../../tasks/skylab/MergeSortBam.wdl" as Merge
import "../../../tasks/skylab/FastqProcessing.wdl" as FastqProcessing
import "../../../tasks/skylab/PairedTagUtils.wdl" as AddBB
import "../../../tasks/broad/Utilities.wdl" as utils
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

    #Expected cells from library preparation
    Int atac_expected_cells = 3000

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

  String pipeline_version = "2.9.3"

  # Determine docker prefix based on cloud provider
  String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
  String acr_docker_prefix = "dsppipelinedev.azurecr.io/"
  String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

  # Docker image names
  String warp_tools_docker = "warp-tools:2.6.1"
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
        input_id = input_id
    }
    if (peak_calling) {
      call peakcalling.PeakCalling as PeakCalling{
        input:
          output_base_name = input_id,
          annotations_gtf = annotations_gtf,
          metrics_h5ad = CreateFragmentFile.Snap_metrics,
          chrom_sizes = chrom_sizes,
          cloud_provider = cloud_provider,
      }
    }
  }
  
  File bam_aligned_output_atac = select_first([BBTag.bb_bam, aligned_bam])
  File fragment_file_atac = select_first([BB_fragment.fragment_file, CreateFragmentFile.fragment_file])
  File fragment_file_index_atac = select_first([BB_fragment.fragment_file_index, CreateFragmentFile.fragment_file_index])
  File snap_metrics_atac = select_first([BB_fragment.Snap_metrics,CreateFragmentFile.Snap_metrics])
  File library_metrics = select_first([BB_fragment.atac_library_metrics, CreateFragmentFile.atac_library_metrics])
    
  output {
    File bam_aligned_output = bam_aligned_output_atac
    File fragment_file = fragment_file_atac
    File fragment_file_index = fragment_file_index_atac
    File snap_metrics = snap_metrics_atac
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
    String gtf_path = annotations_gtf
  }

  parameter_meta {
    bam: "Aligned bam with CB in CB tag. This is the output of the BWAPairedEndAlignment task."
    chrom_sizes: "Text file containing chrom_sizes for genome build (i.e. hg38)."
    annotations_gtf: "GTF for SnapATAC2 to calculate TSS sites of fragment file."
    disk_size: "Disk size used in create fragment file step."
    mem_size: "The size of memory used in create fragment file."
    docker_path: "The docker image path containing the runtime environment for this task"
  }

  command <<<
    set -euo pipefail
    set -x 

    python3 <<CODE

    # import libraries
    import snapatac2.preprocessing as pp
    import snapatac2 as snap
    import scanpy as sc
    import numpy as np
    import polars as pl
    import anndata as ad
    from collections import OrderedDict
    import csv

    # set parameters
    bam = "~{bam}"
    input_id = "~{input_id}"
    chrom_sizes = "~{chrom_sizes}"
    atac_gtf = "~{annotations_gtf}"
    preindex = "~{preindex}"
    atac_nhash_id = "~{atac_nhash_id}"
    mito_list = "~{sep=' ' mito_list}"
    expected_cells = ~{atac_expected_cells}

    print(mito_list)
    mito_list = mito_list.split(" ")
    print("Mitochondrial chromosomes:", mito_list) 
    
    # calculate chrom size dictionary based on text file
    chrom_size_dict={}
    with open('~{chrom_sizes}', 'r') as f:
      for line in f:
        key, value = line.strip().split()
        chrom_size_dict[str(key)] = int(value)

    # extract CB or BB (if preindex is true) tag from bam file to create fragment file
    if preindex == "true":
      data = pp.recipe_10x_metrics("~{bam}", "~{input_id}.fragments.tsv", "temp_metrics.h5ad", is_paired=True, barcode_tag="BB", chrom_sizes=chrom_size_dict, gene_anno=atac_gtf, peaks=None, chrM=mito_list)
    elif preindex == "false":
      data = pp.recipe_10x_metrics("~{bam}", "~{input_id}.fragments.tsv", "temp_metrics.h5ad", is_paired=True, barcode_tag="CB", chrom_sizes=chrom_size_dict, gene_anno=atac_gtf, peaks=None, chrM=mito_list)

    # Add NHashID to metrics 
    data = OrderedDict({'NHashID': atac_nhash_id, **data})
    
    # Calculate atac percent target
    print("Calculating percent target")
    number_of_cells = data['Cells']['Number_of_cells']
    print("Print number of cells", number_of_cells)
    atac_percent_target = number_of_cells / expected_cells*100
    print("Setting percent target in nested dictionary")
    data['Cells']['atac_percent_target'] = atac_percent_target

    # Flatten the dictionary
    flattened_data = []
    for category, metrics in data.items():
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
        writer.writerows(flattened_data)  # Write data

    print(f"Dictionary successfully written to {csv_file_path}")

    atac_data = ad.read_h5ad("temp_metrics.h5ad")
    # Add nhash_id to h5ad file as unstructured metadata
    atac_data.uns['NHashID'] = atac_nhash_id

    # Add GTF to uns field
    # Original path from args.annotation_file
    gtf_path = "~{gtf_path}"  # e.g., 'gs://gcp-public-data--broad-references/hg38/v0/star/v2_7_10a/modified_v43.annotation.gtf'
    
    atac_data.uns["reference_gtf_file"] = gtf_path
    # calculate tsse metrics
    snap.metrics.tsse(atac_data, atac_gtf, exclude_chroms=mito_list)
    # Write new atac file
    atac_data.write_h5ad("~{input_id}.metrics.h5ad")

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
    File atac_library_metrics = "~{input_id}_~{atac_nhash_id}_library_metrics.csv"
  }
}
