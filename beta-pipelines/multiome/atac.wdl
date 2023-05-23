version 1.0

import "../../tasks/skylab/MergeSortBam.wdl" as Merge

workflow ATAC {
  meta {
    description: "Processing for single-cell ATAC-seq data from the level of raw fastq reads. This is the first step of the multiome pipeline. ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a technique used in molecular biology to assess genome-wide chromatin accessibility. This pipeline processes 10x Genomics Multiome ATAC FASTQ files."
    allowNestedInputs: true
  }

  input {
    # Fastq inputs
    Array[File] read1_fastq_gzipped
    Array[File] read2_fastq_gzipped
    Array[File] read3_fastq_gzipped

    # Output prefix/base name for all intermediate files and pipeline outputs
    String output_base_name
    
    # BWA ref 
    File tar_bwa_reference
    
    # CreateFragmentFile input variables 
    Boolean barcodes_in_read_name
    # GTF for SnapATAC2 to calculate TSS sites of fragment file
    File atac_gtf
    # Text file containing chrom_sizes for genome build (i.e. hg38)
    File chrom_sizes

    # script for monitoring tasks 
    File monitoring_script

    # Whitelist
    File whitelist

    # TrimAdapters input 
    String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
  }

  parameter_meta {
    read1_fastq_gzipped: "read 1 FASTQ file as input for the pipeline, contains read 1 of paired reads"
    read2_fastq_gzipped: "read 2 FASTQ file as input for the pipeline, contains the cellular barcodes corresponding to the reads in the read1 FASTQ and read 3 FASTQ"
    read3_fastq_gzipped: "read 3 FASTQ file as input for the pipeline, contains read 2 of paired reads"
    output_base_name: "base name to be used for the pipelines output and intermediate files"
    monitoring_script : "script to monitor resource comsumption of tasks"
    tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"

  }

  call FastqProcessing as SplitFastq {
    input:
      read1_fastq = read1_fastq_gzipped,
      read3_fastq = read3_fastq_gzipped,
      barcodes_fastq = read2_fastq_gzipped,
      output_base_name = output_base_name,
      whitelist = whitelist
   }

  scatter(idx in range(length(SplitFastq.fastq_R1_output_array))) {

    call TrimAdapters {
      input:
        read1_fastq = SplitFastq.fastq_R1_output_array[idx], 
        read3_fastq = SplitFastq.fastq_R3_output_array[idx],
        output_base_name = output_base_name + "_" + idx,
        monitoring_script = monitoring_script,
        adapter_seq_read1 = adapter_seq_read1,
        adapter_seq_read3 = adapter_seq_read3
    }

    call BWAPairedEndAlignment {
      input:
        read1_fastq = TrimAdapters.fastq_trimmed_adapter_output_read1,
        read3_fastq = TrimAdapters.fastq_trimmed_adapter_output_read3,
        tar_bwa_reference = tar_bwa_reference,
        output_base_name = output_base_name + "_" + idx,
        monitoring_script = monitoring_script
      }
      
    call AddCBtags {
     input:
       bam = BWAPairedEndAlignment.bam_aligned_output,
       output_base_name = output_base_name
   }  
  
  } 

    call Merge.MergeSortBamFiles as MergeBam {
    input:
      bam_inputs = AddCBtags.output_cb_bam,
      output_bam_filename = output_base_name + ".bam",
      sort_order = "coordinate"
  }
  

  call CreateFragmentFile {
    input:
      bam = MergeBam.output_bam,
      chrom_sizes = chrom_sizes,
      barcodes_in_read_name = barcodes_in_read_name,
      atac_gtf = atac_gtf
   }

  output {
    File bam_aligned_output = MergeBam.output_bam
    File fragment_file = CreateFragmentFile.fragment_file
    File snap_metrics = CreateFragmentFile.Snap_metrics
  }
}

task FastqProcessing {

  input {
    Array[File] read1_fastq
    Array[File] read3_fastq
    Array[File] barcodes_fastq
    String read_structure = "16C"
    String barcode_orientation = "FIRST_BP_RC"
    String output_base_name
    File whitelist

    # [?] copied from corresponding optimus wdl for fastqprocessing 
    # using the latest build of warp-tools in GCR
    String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1682971351"
    # Runtime attributes [?]
    Int mem_size = 5
    Int cpu = 16   
    # TODO decided cpu
    # estimate that bam is approximately equal in size to fastq, add 20% buffer
    Int disk_size = ceil(2 * ( size(read1_fastq, "GiB") + size(read3_fastq, "GiB") + size(barcodes_fastq, "GiB") )) + 400
    Int preemptible = 3
   }

  meta {
    description: "Converts a set of fastq files to unaligned bam file/fastq file, also corrects barcodes and partitions the alignments by barcodes. Allows for variable barcode and umi lengths as input, if applicable."
   }

  parameter_meta {
        read1_fastq: "Array of read 1 FASTQ files of paired reads -- forward reads"
        read3_fastq: "Array of read 3 FASTQ files of paired reads -- reverse reads"
        barcodes_fastq: "Array of read 2 FASTQ files which contains the cellular barcodes"
        output_base_name: "Name of sample matching this file, inserted into read group header"
        read_structure: "A string that specifies the barcode (C) positions in the Read 2 fastq"
        barcode_orientation: "A string that specifies the orientation of barcode needed for scATAC data. The default is FIRST_BP. Other options include LAST_BP, FIRST_BP_RC or LAST_BP_RC."
        whitelist: "10x genomics cell barcode whitelist for scATAC"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_size: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        disk_size: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
   }

  command <<<
  
    set -euo pipefail
    
    # Cat files for each r1, r3 and barcodes together 
    cat ~{sep=' ' read1_fastq} > r1.fastq.gz
    cat ~{sep=' ' read3_fastq} > r3.fastq.gz
    cat ~{sep=' ' barcodes_fastq} > barcodes.fastq.gz
    
    # Call fastq process 
    # outputs fastq files where the corrected barcode is in the read name
    fastqprocess \
        --bam-size 10.0 \
        --sample-id "~{output_base_name}" \
        --R1 barcodes.fastq.gz \
        --R2 r1.fastq.gz \
        --R3 r3.fastq.gz \
        --white-list "~{whitelist}" \
        --output-format "FASTQ" \
        --barcode-orientation "~{barcode_orientation}" \
        --read-structure "~{read_structure}"
        
    >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem_size} MiB"
    disks: "local-disk ${disk_size} HDD"
    preemptible: preemptible
   }
  
  output {
    Array[File] fastq_R1_output_array = glob("fastq_R1_*")
    Array[File] fastq_R3_output_array = glob("fastq_R3_*")
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
      String docker_image = "quay.io/broadinstitute/cutadapt:1.18"
      File monitoring_script
   }

   parameter_meta {
      read1_fastq: "read 1 fastq file containing sequencing reads as input for the pipeline"
      read3_fastq: "read 3 fastq file containing sequencing reads as input for the pipeline"
      min_length: "the minimum length for trimming. Reads that are too short even before adapter removal are also discarded"
      quality_cutoff: "cutadapt option to trim low-quality ends from reads before adapter removal"
      adapter_seq_read1: "cutadapt option for the sequence adapter for read 1 fastq"
      adapter_seq_read3: "cutadapt option for the sequence adapter for read 3 fastq"
      output_base_name: "base name to be used for the output of the task"
      docker_image: "the docker image using cutadapt to be used (default: quay.io/broadinstitute/cutadapt:1.18)"
      monitoring_script : "script to monitor resource consumption of tasks"
      mem_size: "the size of memory used during trimming adapters"
      disk_size : "disk size used in trimming adapters step"
   }
      
    # output names for trimmed reads
    String fastq_trimmed_adapter_output_name_read1 = output_base_name + ".R1.trimmed_adapters.fastq.gz"
    String fastq_trimmed_adapter_output_name_read3 = output_base_name + ".R3.trimmed_adapters.fastq.gz"
   
    # using cutadapt to trim off sequence adapters
    command <<<
      set -euo pipefail

      if [ ! -z "~{monitoring_script}" ]; then
        chmod a+x ~{monitoring_script}
        ~{monitoring_script} > monitoring.log &
      else
        echo "No monitoring script given as input" > monitoring.log &
      fi

      # fastq's, "-f", -A for paired adapters read 2"
      cutadapt \
        -f fastq \
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
      File monitoring_log = "monitoring.log"
    }
  }

# align the two trimmed fastq as paired end data using BWA
task BWAPairedEndAlignment {
    input {
      File read1_fastq
      File read3_fastq
      File tar_bwa_reference
      String read_group_id = "RG1"
      String read_group_sample_name = "RGSN1"
      String output_base_name
      String docker_image = "us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091"
      
      # script for monitoring tasks 
      File monitoring_script
      
      # Runtime attributes
      Int disk_size = ceil(3.25 * (size(read1_fastq, "GiB") + size(read3_fastq, "GiB") + size(tar_bwa_reference, "GiB"))) + 200 
      Int nthreads = 16
      Int mem_size = 8
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
      docker_image: "the docker image using BWA to be used (default: us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091)"
      monitoring_script : "script to monitor resource comsumption of tasks"
    }

    String bam_aligned_output_name = output_base_name + ".aligned.bam"

    # bwa and call samtools to convert sam to bam
    command <<<

      set -euo pipefail

      if [ ! -z "~{monitoring_script}" ]; then
        chmod a+x ~{monitoring_script}
         ~{monitoring_script} > monitoring.log &
      else
        echo "No monitoring script given as input" > monitoring.log &
      fi

      # prepare reference
      declare -r REF_DIR=$(mktemp -d genome_referenceXXXXXX)
      tar -xf "~{tar_bwa_reference}" -C $REF_DIR --strip-components 1
      rm "~{tar_bwa_reference}"

      # align w/ BWA: -t for number of cores
      bwa \
        mem \
        -R "@RG\tID:~{read_group_id}\tSM:~{read_group_sample_name}" \
        -t ~{nthreads} \
        $REF_DIR/genome.fa \
        ~{read1_fastq} ~{read3_fastq} \
        | samtools view -bS - > ~{bam_aligned_output_name}    
     >>>

    runtime {
      docker: docker_image
      disks: "local-disk ${disk_size} HDD"
      cpu: nthreads
      memory: "${mem_size} GiB"
    }

    output {
      File bam_aligned_output = bam_aligned_output_name
      File monitoring_log = "monitoring.log"
    }
  }

#add CB and CR tags to BAM file
task AddCBtags {
  input {
    File bam
    String output_base_name
    Int disk_size = ceil(size(bam, "GiB") + 50)
    Int mem_size = 10
  }
  
  String bam_cb_output_name = output_base_name + ".cb.aligned.bam"
  
  command <<<
    # We reused tags from Broad's Share-seq pipeline; XC tag is not currently needed by pipeline 
    samtools view -h ~{bam} \
    | awk '{if ($0 ~ /^@/) {print $0} else {cr=substr($1, index($1, "CR:")+3); n=split($1, a,":CB:"); if (n == 2) {cb="CB:Z:"a[1]"\t";} else {cb="";} print($0 "\tCR:Z:" cr "\t" cb "XC:Z:" substr($1, index($1, "CB:")+3));}}' \
    | samtools view -b -o ~{bam_cb_output_name}
    
    # Piping to samtools sort works, but isn't necessary for SnapATAC2
    #| \ samtools sort -o ~{bam_cb_output_name}
  >>>

  output {
    File output_cb_bam = "~{bam_cb_output_name}"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091"
    disks: "local-disk ${disk_size} HDD"
    memory: "${mem_size} GiB"
  }
}
# make fragment file
task CreateFragmentFile {
  input {
    File bam
    File atac_gtf
    File chrom_sizes
    Boolean barcodes_in_read_name
    Int disk_size = ceil(size(bam, "GiB") + 200)
    Int mem_size = 50
  }

  String bam_base_name = basename(bam, ".bam")

  parameter_meta {
    bam: "Aligned bam with CB in CB tag. This is the output of the AddCBtags task."
    atac_gtf: "GTF for SnapATAC2 to calculate TSS sites of fragment file."
    chrom_sizes: "Text file containing chrom_sizes for genome build (i.e. hg38)."
    barcodes_in_read_name: "Set to True if barcodes are in the read names, and False otherwise. "
    disk_size: "Disk size used in create fragment file step."
    mem_size: "The size of memory used in create fragment file."
  }

  command <<<
    set -e pipefail

    python3 <<CODE

    # set parameters
    atac_gtf = "~{atac_gtf}"
    barcodes_in_read_name = "~{barcodes_in_read_name}"
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

    # if barcodes are in the read name, then use barcode_regex to extract them. otherwise, use barcode_tag
    if barcodes_in_read_name=="true":
      pp.make_fragment_file("~{bam}", "~{bam_base_name}.fragments.tsv", is_paired=True, barcode_regex="([^:]*)")
    elif barcodes_in_read_name=="false":
      pp.make_fragment_file("~{bam}", "~{bam_base_name}.fragments.tsv", is_paired=True, barcode_tag="CB")
    
    # calculate quality metrics; note min_num_fragments and min_tsse are set to 0 instead of default
    # those settings allow us to retain all barcodes
    pp.import_data("~{bam_base_name}.fragments.tsv", file="~{bam_base_name}.metrics.h5ad", chrom_size=chrom_size_dict, gene_anno="~{atac_gtf}", min_num_fragments=0, min_tsse=0)

    CODE
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/snapatac2:1.0.3-2.3.0-1682089891"
    disks: "local-disk ${disk_size} HDD"
    memory: "${mem_size} GiB"
  }

  output {
    File fragment_file = "~{bam_base_name}.fragments.tsv"
    File Snap_metrics = "~{bam_base_name}.metrics.h5ad"
  }
}



