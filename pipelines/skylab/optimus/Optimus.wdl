version 1.0

import "../../../tasks/skylab/FastqProcessing.wdl" as FastqProcessing
import "../../../tasks/skylab/StarAlign.wdl" as StarAlign
import "../../../tasks/skylab/Metrics.wdl" as Metrics
import "../../../tasks/skylab/RunEmptyDrops.wdl" as RunEmptyDrops
import "../../../tasks/skylab/LoomUtils.wdl" as LoomUtils
import "../../../tasks/skylab/CheckInputs.wdl" as OptimusInputChecks
import "../../../tasks/skylab/MergeSortBam.wdl" as Merge

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

    # Chemistry options include: 2 or 3
    Int tenx_chemistry_version
    # Whitelist is selected based on the input tenx_chemistry_version
    File whitelist = checkOptimusInput.whitelist_out

    # Emptydrops lower cutoff
    Int emptydrops_lower = 100

    # Set to true to override input checks and allow pipeline to proceed with invalid input
    Boolean force_no_check = false
    
    # Check that tenx_chemistry_version matches the length of the read 1 fastq;
    # Set to true if you expect that r1_read_length does not match length of UMIs/barcodes for 10x chemistry v2 (26 bp) or v3 (28 bp).
    Boolean ignore_r1_read_length = false

    # Set to true to count reads in stranded mode
    String use_strand_info = "false"
    
# Set to true to count reads aligned to exonic regions in sn_rna mode
    Boolean count_exons = false

    # this pipeline does not set any preemptible varibles and only relies on the task-level preemptible settings
    # you could override the tasklevel preemptible settings by passing it as one of the workflows inputs
    # for example: `"Optimus.StarAlign.preemptible": 3` will let the StarAlign task, which by default disables the
    # usage of preemptible machines, attempt to request for preemptible instance up to 3 times. 
  }

  # version of this pipeline
  
  String pipeline_version = "5.7.0"

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
    use_strand_info: "Set to true to count reads in stranded mode"
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

  call FastqProcessing.FastqProcessing as SplitFastq {
    input:
      i1_fastq = i1_fastq,
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      whitelist = whitelist,
      chemistry = tenx_chemistry_version,
      sample_id = input_id
  }

  scatter(idx in range(length(SplitFastq.fastq_R1_output_array))) {
    call StarAlign.STARsoloFastq as STARsoloFastq {
      input:
        r1_fastq = [SplitFastq.fastq_R1_output_array[idx]],
        r2_fastq = [SplitFastq.fastq_R2_output_array[idx]],
        white_list = whitelist,
        tar_star_reference = tar_star_reference,
        chemistry = tenx_chemistry_version,
        counting_mode = counting_mode,
        count_exons = count_exons,
        output_bam_basename = output_bam_basename + "_" + idx
    }
  }
  call Merge.MergeSortBamFiles as MergeBam {
    input:
      bam_inputs = STARsoloFastq.bam_output,
      output_bam_filename = output_bam_basename + ".bam",
      sort_order = "coordinate"
  }
  call Metrics.CalculateGeneMetrics as GeneMetrics {
    input:
      bam_input = MergeBam.output_bam,
      mt_genes = mt_genes,
      input_id = input_id
  }

  call Metrics.CalculateCellMetrics as CellMetrics {
    input:
      bam_input = MergeBam.output_bam,
      mt_genes = mt_genes,
      original_gtf = annotations_gtf,
      input_id = input_id
  }

  call StarAlign.MergeStarOutput as MergeStarOutputs {
    input:
      barcodes = STARsoloFastq.barcodes,
      features = STARsoloFastq.features,
      matrix = STARsoloFastq.matrix,
      input_id = input_id
  }
  if (counting_mode == "sc_rna"){
    call RunEmptyDrops.RunEmptyDrops {
      input:
        sparse_count_matrix = MergeStarOutputs.sparse_counts,
        row_index = MergeStarOutputs.row_index,
        col_index = MergeStarOutputs.col_index,
        emptydrops_lower = emptydrops_lower
    }
  }

  if (!count_exons) {
    call LoomUtils.OptimusLoomGeneration{
      input:
        input_id = input_id,
        input_name = input_name,
        input_id_metadata_field = input_id_metadata_field,
        input_name_metadata_field = input_name_metadata_field,
        annotation_file = annotations_gtf,
        cell_metrics = CellMetrics.cell_metrics,
        gene_metrics = GeneMetrics.gene_metrics,
        sparse_count_matrix = MergeStarOutputs.sparse_counts,
        cell_id = MergeStarOutputs.row_index,
        gene_id = MergeStarOutputs.col_index,
        empty_drops_result = RunEmptyDrops.empty_drops_result,
        counting_mode = counting_mode,
        pipeline_version = "Optimus_v~{pipeline_version}"
    }
  }
  if (count_exons  && counting_mode=="sn_rna") {
    call StarAlign.MergeStarOutput as MergeStarOutputsExons {
      input:
        barcodes = STARsoloFastq.barcodes_sn_rna,
        features = STARsoloFastq.features_sn_rna,
        matrix = STARsoloFastq.matrix_sn_rna,
        input_id = input_id
    }
    call LoomUtils.SingleNucleusOptimusLoomOutput as OptimusLoomGenerationWithExons{
      input:
        input_id = input_id,
        input_name = input_name,
        input_id_metadata_field = input_id_metadata_field,
        input_name_metadata_field = input_name_metadata_field,
        annotation_file = annotations_gtf,
        cell_metrics = CellMetrics.cell_metrics,
        gene_metrics = GeneMetrics.gene_metrics,
        sparse_count_matrix = MergeStarOutputs.sparse_counts,
        cell_id = MergeStarOutputs.row_index,
        gene_id = MergeStarOutputs.col_index,
        sparse_count_matrix_exon = MergeStarOutputsExons.sparse_counts,
        cell_id_exon = MergeStarOutputsExons.row_index,
        gene_id_exon = MergeStarOutputsExons.col_index,
        pipeline_version = "Optimus_v~{pipeline_version}"
    }

  }

  File final_loom_output = select_first([OptimusLoomGenerationWithExons.loom_output, OptimusLoomGeneration.loom_output])


  output {
    # version of this pipeline
    String pipeline_version_out = pipeline_version
    File bam = MergeBam.output_bam
    File matrix = MergeStarOutputs.sparse_counts
    File matrix_row_index = MergeStarOutputs.row_index
    File matrix_col_index = MergeStarOutputs.col_index
    File cell_metrics = CellMetrics.cell_metrics
    File gene_metrics = GeneMetrics.gene_metrics
    File? cell_calls = RunEmptyDrops.empty_drops_result
    # loom
    File loom_output_file = final_loom_output
}
}
