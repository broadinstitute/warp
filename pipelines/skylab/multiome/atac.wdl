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

    # Whitelist
    File whitelist

    # TrimAdapters input
    String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
  }

  String pipeline_version = "1.0.0"

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
      whitelist = whitelist
  }

  call Bowtie2Build {
    input:
      fasta_input = "gs://gcp-public-data--broad-references/hg38/v0/GRCh38.primary_assembly.genome.fa",
      index_prefix = "GRCh38"
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

    call share_atac_align {
      input:
        fastq_R1 = TrimAdapters.fastq_trimmed_adapter_output_read1,
        fastq_R2 = TrimAdapters.fastq_trimmed_adapter_output_read3,
        multimappers = 5,
        genome_index_tar = Bowtie2Build.bowtie2_index_files,
        chemistry = "atac",
        genome_name = "GRCh39",
        prefix = "bowtie2test",
    }
  }

  call Merge.MergeSortBamFiles as MergeBam {
    input:
      output_bam_filename = input_id + ".bam",
      bam_inputs = share_atac_align.atac_alignment,
      sort_order = "coordinate"
  }


  call CreateFragmentFile {
    input:
      bam = MergeBam.output_bam,
      chrom_sizes = chrom_sizes,
      annotations_gtf = annotations_gtf
  }

  output {
    File bam_aligned_output = MergeBam.output_bam
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

## align the two trimmed fastq as paired end data using BWA
#task BWAPairedEndAlignment {
#  input {
#    File read1_fastq
#    File read3_fastq
#    File tar_bwa_reference
#    String read_group_id = "RG1"
#    String read_group_sample_name = "RGSN1"
#    String output_base_name
#    String docker_image = "us.gcr.io/broad-gotc-prod/samtools-bwa-mem-2:1.0.0-2.2.1_x64-linux-1685469504"
#
#    # Runtime attributes
#    Int disk_size = ceil(3.25 * (size(read1_fastq, "GiB") + size(read3_fastq, "GiB") + size(tar_bwa_reference, "GiB"))) + 400
#    Int nthreads = 16
#    Int mem_size = 40
#  }
#
#  parameter_meta {
#    read1_fastq: "the trimmed read 1 fastq file containing sequencing reads as input for the aligner"
#    read3_fastq: "the trimmed read 1 fastq file containing sequencing reads as input for the aligner"
#    tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"
#    read_group_id: "the read group id to be added upon alignment"
#    read_group_sample_name: "the read group sample to be added upon alignment"
#    nthreads: "the number of threads to use during bwa alignment"
#    mem_size: "the size of memory used during alignment"
#    disk_size : "disk size used in bwa alignment step"
#    output_base_name: "basename to be used for the output of the task"
#    docker_image: "the docker image using BWA to be used (default: us.gcr.io/broad-gotc-prod/samtools-bwa-mem-2:1.0.0-2.2.1_x64-linux-1685469504)"
#  }
#
#  String bam_aligned_output_name = output_base_name + ".aligned.bam"
#
#  # bwa and call samtools to convert sam to bam
#  command <<<
#
#    set -euo pipefail
#
#    # prepare reference
#    declare -r REF_DIR=$(mktemp -d genome_referenceXXXXXX)
#    tar -xf "~{tar_bwa_reference}" -C $REF_DIR --strip-components 1
#    rm "~{tar_bwa_reference}"
#
#    # align w/ BWA: -t for number of cores
#    bwa-mem2 \
#    mem \
#    -R "@RG\tID:~{read_group_id}\tSM:~{read_group_sample_name}" \
#    -C \
#    -t ~{nthreads} \
#    $REF_DIR/genome.fa \
#    ~{read1_fastq} ~{read3_fastq} \
#    | samtools view -bS - > ~{bam_aligned_output_name}
#  >>>
#
#  runtime {
#    docker: docker_image
#    disks: "local-disk ${disk_size} HDD"
#    cpu: nthreads
#    memory: "${mem_size} GiB"
#  }
#
#  output {
#    File bam_aligned_output = bam_aligned_output_name
#  }
#}
# align the two trimmed fastq as paired end data using BWA

# make fragment file
task CreateFragmentFile {
  input {
    File bam
    File annotations_gtf
    File chrom_sizes
    Int disk_size = ceil(size(bam, "GiB") + 200)
    Int mem_size = 50
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
    disks: "local-disk ${disk_size} HDD"
    memory: "${mem_size} GiB"
  }

  output {
    File fragment_file = "~{bam_base_name}.fragments.tsv"
    File Snap_metrics = "~{bam_base_name}.metrics.h5ad"
  }
}


# build with the converted fasta's
task Bowtie2Build {
  input {
    File fasta_input
    String index_prefix
  }

  # input file size
  Float input_size = size(fasta_input, "GB")

  # use bowtie2 to build indexes
  command <<<

    bowtie2-build \
    -f ~{fasta_input} ~{index_prefix}

    tar -zcvf ~{index_prefix}.index.tar.gz ~{index_prefix}*

  >>>

  runtime {
    docker: "quay.io/broadinstitute/bowtie2:2.3.4.3"
    disks: "local-disk 50 HDD"
    cpu: 1
    memory: "20 GB"
  }

  output {
    File bowtie2_index_files = "~{index_prefix}.index.tar.gz"
  }
}


task share_atac_align {
  meta {
    version: 'v0.1'
    author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
    description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: align ATAC task using bowtie2'
  }

  input {
    # This task takes in input the preprocessed ATAC fastqs and align them to the genome.
    File fastq_R1
    File fastq_R2
    Int? multimappers # = 5
    File genome_index_tar       # This is a tar.gz folder with all the index files.
    String chemistry
    String genome_name          # GRCh38, mm10
    String prefix = "sample-share"
    Int? cpus = 16
    Float? disk_factor = 8.0
    Float? memory_factor = 0.15
    String? docker_image = "us.gcr.io/buenrostro-share-seq/share_task_bowtie2"
  }

  # Determine the size of the input
  Float input_file_size_gb = size(fastq_R1, "G") + size(fastq_R2, "G")

  # Determining memory size base on the size of the input files.
  Float mem_gb = 5.0 + size(genome_index_tar, "G") + memory_factor * input_file_size_gb

  # Determining disk size base on the size of the input files.
  Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

  # Determining disk type base on the size of disk.
  String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

  # Determining memory for samtools.
  Float samtools_memory_gb = 0.8 * mem_gb # Samtools has overheads so reducing the memory to 80% of the total.

  # Number of threads to beable to use 4GB of memory per thread seems to be the fastest way
  Int samtools_threads_ = floor(samtools_memory_gb / 4)
  Int samtools_threads =  if samtools_threads_ == 0 then 1 else samtools_threads_

  # Now that we know how many threads we can use to assure 4GB of memory per thread
  # we assign any remaining memory to the threads.
  Int samtools_memory_per_thread_ = floor(samtools_memory_gb * 1024 / samtools_threads) # Computing the memory per thread for samtools in MB.
  Int samtools_memory_per_thread = if samtools_memory_per_thread_ < 768 then 768 else samtools_memory_per_thread_


  # Define tmp file name
  String unsorted_bam = "${prefix}.atac.align.${genome_name}.bam"

  # Define the output names
  String sorted_bam = "${prefix}.atac.align.k${multimappers}.${genome_name}.sorted.bam"
  String sorted_bai = "${prefix}.atac.align.k${multimappers}.${genome_name}.sorted.bam.bai"
  String alignment_log = "${prefix}.atac.align.k${multimappers}.${genome_name}.log"


  command <<<
    set -e

    tar zxvf ~{genome_index_tar} --no-same-owner -C ./
    genome_prefix=$(basename $(find . -type f -name "*.rev.1.bt2") .rev.1.bt2)

    # Aligning and adding the cell barcode to the CB tag and the barcodes plus pkr in the XC tag.
    bowtie2 -X2000 \
        -p ~{cpus} \
        ~{if defined(multimappers) then "-k ~{multimappers}" else ""} --rg-id ~{prefix + "."}atac \
        --rg "SM:None" \
        --rg "LB:None" \
        --rg "PL:Illumina" \
        ~{if "~{chemistry}" != "shareseq" then "--sam-append-comment" else ""} \
        -x $genome_prefix \
        -1 ~{fastq_R1} \
        -2 ~{fastq_R2} 2> ~{alignment_log} | \
        samtools view \
            -b \
            -S \
            -@ ~{samtools_threads} \
            - \
            -o ~{unsorted_bam}


    #if [ '~{chemistry}' != 'shareseq' ]; then
    #    samtools sort \
    #        -@ ~{samtools_threads} \
    #        -m ~{samtools_memory_per_thread}M \
    #        ~{unsorted_bam} \
    #        -o ~{sorted_bam}
    #else
    #    # Splitting the read name to ge the cell barcode and adding it to the CB tag in the BAM file.
    #    samtools view -h ~{unsorted_bam} | \
    #    awk '{if ($0 ~ /^@/) {print $0} else {split($1,a,"[,_]"); print($0 "\tCB:Z:" a[2]a[3]a[4] "\tXC:Z:" a[2]a[3]a[4] "_" a[5]);}}' | \
    #    samtools sort \
    #        -@ ~{samtools_threads} \
    #        -m ~{samtools_memory_per_thread}M \
    #        - \
    #        -o ~{sorted_bam}
    #fi

    #samtools index -@ ~{cpus} ~{sorted_bam}

  >>>

  output {
    File atac_alignment = unsorted_bam
    #File atac_alignment_index = sorted_bai
    File atac_alignment_log = alignment_log
  }

  runtime {
    cpu: cpus
    docker: "${docker_image}"
    disks: "local-disk ${disk_gb} ${disk_type}"
    maxRetries:1
    memory: "${mem_gb} GB"
    memory_retry_multiplier: 2
  }
}