version 1.0

import "../../../tasks/skylab/FastqProcessing.wdl" as FastqProcessing
import "../../../tasks/skylab/StarAlign.wdl" as StarAlign
import "../../../tasks/skylab/Metrics.wdl" as Metrics
import "../../../tasks/skylab/RunEmptyDrops.wdl" as RunEmptyDrops
import "../../../tasks/skylab/CheckInputs.wdl" as OptimusInputChecks
import "../../../tasks/skylab/MergeSortBam.wdl" as Merge
import "../../../tasks/skylab/H5adUtils.wdl" as H5adUtils
import "../../../tasks/skylab/GetSplits.wdl" as GetSplits

workflow Optimus {
  meta {
    description: "The optimus 3' pipeline processes 10x genomics sequencing data based on the v2 chemistry. It corrects cell barcodes and UMIs, aligns reads, marks duplicates, and returns data as alignments in BAM format and as counts in sparse matrix exchange format."
  }

  input {
    # Mode for counting either "sc_rna" or "sn_rna"
    String counting_mode = "sc_rna"

    # Sequencing data inputs
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    String input_id
    String output_bam_basename = input_id
    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field
    # organism reference parameters
    File tar_star_reference
    File annotations_gtf
    File ref_genome_fasta
    File? mt_genes
    String? soloMultiMappers

    # Chemistry options include: 2 or 3
    Int tenx_chemistry_version
    # Whitelist is selected based on the input tenx_chemistry_version
    File whitelist = checkOptimusInput.whitelist_out

    # read_structure is based on v2 or v3 chemistry
    String read_struct = checkOptimusInput.read_struct_out

    # Emptydrops lower cutoff
    Int emptydrops_lower = 100

    # Set to true to override input checks and allow pipeline to proceed with invalid input
    Boolean force_no_check = false
    
    # Check that tenx_chemistry_version matches the length of the read 1 fastq;
    # Set to true if you expect that r1_read_length does not match length of UMIs/barcodes for 10x chemistry v2 (26 bp) or v3 (28 bp).
    Boolean ignore_r1_read_length = false

    # Set to Forward, Reverse, or Unstranded to account for stranded library preparations (per STARsolo documentation)
    String star_strand_mode = "Forward"
    
    # Set to true to count reads aligned to exonic regions in sn_rna mode
    Boolean count_exons = false
     

    # Star machine type -- to select number of splits 
    Int num_threads_star = 128
    Int mem_size_star = 512
    String cpu_platform_star = "Intel Ice Lake"
 
    # this pipeline does not set any preemptible varibles and only relies on the task-level preemptible settings
    # you could override the tasklevel preemptible settings by passing it as one of the workflows inputs
    # for example: `"Optimus.StarAlign.preemptible": 3` will let the StarAlign task, which by default disables the
    # usage of preemptible machines, attempt to request for preemptible instance up to 3 times. 
  }

  # version of this pipeline

  String pipeline_version = "6.3.5"

  # this is used to scatter matched [r1_fastq, r2_fastq, i1_fastq] arrays
  Array[Int] indices = range(length(r1_fastq))

  # 10x parameters
  File whitelist_v2 = "gs://gcp-public-data--broad-references/RNA/resources/737k-august-2016.txt"
  File whitelist_v3 = "gs://gcp-public-data--broad-references/RNA/resources/3M-febrary-2018.txt"
  # Takes the first read1 FASTQ from the inputs to check for chemistry match
  File r1_single_fastq = r1_fastq[0]

  parameter_meta {
    r1_fastq: "forward read, contains cell barcodes and molecule barcodes"
    r2_fastq: "reverse read, contains cDNA fragment generated from captured mRNA"
    i1_fastq: "(optional) index read, for demultiplexing of multiple samples on one flow cell."
    input_id: "name of sample matching this file, inserted into read group header"
    input_id_metadata_field: "String that describes the metadata field containing the input_id"
    input_name: "User provided sample name or cell_names"
    input_name_metadata_field: "String that describes the metadata field containing the input_name"
    tar_star_reference: "star genome reference"
    annotations_gtf: "gtf containing annotations for gene tagging (must match star reference)"
    ref_genome_fasta: "genome fasta file (must match star reference)"
    whitelist: "10x genomics cell barcode whitelist"
    tenx_chemistry_version: "10X Genomics v2 (10 bp UMI) or v3 chemistry (12bp UMI)"
    force_no_check: "Set to true to override input checks and allow pipeline to proceed with invalid input"
    star_strand_mode: "STAR mode for handling stranded reads. Options are 'Forward', 'Reverse, or 'Unstranded.' Default is Forward."
  }

  call OptimusInputChecks.checkOptimusInput {
    input:
      force_no_check = force_no_check,
      counting_mode = counting_mode,
      count_exons = count_exons,
      whitelist_v2 = whitelist_v2,
      whitelist_v3 = whitelist_v3,
      tenx_chemistry_version = tenx_chemistry_version,
      r1_fastq = r1_single_fastq,
      ignore_r1_read_length = ignore_r1_read_length
  }

  call StarAlign.STARGenomeRefVersion as ReferenceCheck {
    input:
      tar_star_reference = tar_star_reference
  }
  
  call GetSplits.GetNumSplits as GetNumSplits {
    input:
       nthreads = num_threads_star, 
       mem_size = mem_size_star,
       cpu_platform = cpu_platform_star
  }

  call FastqProcessing.FastqProcessing as SplitFastq {
    input:
      i1_fastq = i1_fastq,
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      whitelist = whitelist,
      chemistry = tenx_chemistry_version,
      sample_id = input_id,
      read_struct = read_struct,
      num_output_files = GetNumSplits.ranks_per_node_out
  }

  call StarAlign.STARsoloFastqTest as STARsoloFastq {
      input:
        r1_fastq = SplitFastq.fastq_R1_output_array,
        r2_fastq = SplitFastq.fastq_R2_output_array,
        star_strand_mode = star_strand_mode,
        white_list = whitelist,
        tar_star_reference = tar_star_reference,
        chemistry = tenx_chemistry_version,
        counting_mode = counting_mode,
        count_exons = count_exons,
        output_bam_basename = output_bam_basename,
        soloMultiMappers = soloMultiMappers,
        nthreads = num_threads_star, 
        mem_size = mem_size_star,
        cpu_platform = cpu_platform_star
    }
  
  # call Merge.MergeSortBamFiles as MergeBam {
  #   input:
  #     bam_inputs = STARsoloFastq.bam_aligned_output,
  #     output_bam_filename = output_bam_basename + ".bam",
  #     sort_order = "coordinate"
  # }
  call Metrics.CalculateGeneMetrics as GeneMetrics {
    input:
      bam_input = STARsoloFastq.bam_aligned_output,
      mt_genes = mt_genes,
      input_id = input_id
  }

  call Metrics.CalculateCellMetrics as CellMetrics {
    input:
      bam_input = STARsoloFastq.bam_aligned_output,
      mt_genes = mt_genes,
      original_gtf = annotations_gtf,
      input_id = input_id
  }

}
