version 1.0

import "FastqToUBam.wdl" as FastqToUBam
import "Attach10xBarcodes.wdl" as Attach
import "SplitBamByCellBarcode.wdl" as Split
import "MergeSortBam.wdl" as Merge
import "CreateCountMatrix.wdl" as Count
import "StarAlignBamSingleEnd.wdl" as StarAlignBam
import "TagGeneExon.wdl" as TagGeneExon
import "SequenceDataWithMoleculeTagMetrics.wdl" as Metrics
import "TagSortBam.wdl" as TagSortBam
import "RunEmptyDrops.wdl" as RunEmptyDrops
import "LoomUtils.wdl" as LoomUtils
import "Picard.wdl" as Picard
import "UmiCorrection.wdl" as UmiCorrection
import "ScatterBam.wdl" as ScatterBam
import "ModifyGtf.wdl" as ModifyGtf
import "OptimusInputChecks.wdl" as OptimusInputChecks

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
    String sample_id
    String? output_bam_basename = sample_id

    # organism reference parameters
    File tar_star_reference
    File annotations_gtf
    File ref_genome_fasta

    # 10x parameters
    File whitelist
    # tenX_v2, tenX_v3
    String chemistry = "tenX_v2" 

    # environment-specific parameters
    String fastq_suffix = ""
    
    # Emptydrops lower cutoff
    Int emptydrops_lower = 100

    # Set to true to override input checks and allow pipeline to proceed with invalid input
    Boolean force_no_check = false

    # this pipeline does not set any preemptible varibles and only relies on the task-level preemptible settings
    # you could override the tasklevel preemptible settings by passing it as one of the workflows inputs
    # for example: `"Optimus.StarAlign.preemptible": 3` will let the StarAlign task, which by default disables the
    # usage of preemptible machines, attempt to request for preemptible instance up to 3 times. 
  }

  # version of this pipeline
  String version = "optimus_v3.0.0"

  # this is used to scatter matched [r1_fastq, r2_fastq, i1_fastq] arrays
  Array[Int] indices = range(length(r1_fastq))

  parameter_meta {
    r1_fastq: "forward read, contains cell barcodes and molecule barcodes"
    r2_fastq: "reverse read, contains cDNA fragment generated from captured mRNA"
    i1_fastq: "(optional) index read, for demultiplexing of multiple samples on one flow cell."
    sample_id: "name of sample matching this file, inserted into read group header"
    tar_star_reference: "star genome reference"
    annotations_gtf: "gtf containing annotations for gene tagging (must match star reference)"
    ref_genome_fasta: "genome fasta file (must match star reference)"
    whitelist: "10x genomics cell barcode whitelist"
    tenX_v3_chemistry: "assume 10X Genomics v3 chemistry with 12bp UMI (in contrast to default v2 with 10bp UMI)"
    fastq_suffix: "when running in green box, need to add '.gz' for picard to detect the compression"
    output_zarr: "whether to run the taks that converts the outputs to Zarr format, by default it's true"
    force_no_check: "Set to true to override input checks and allow pipeline to proceed with invalid input"
  }

  call OptimusInputChecks.checkOptimusInput {
    input:
      force_no_check = force_no_check,
      chemistry = chemistry,
      counting_mode = counting_mode
  }

  scatter (index in indices) {
    call FastqToUBam.FastqToUBam {
      input:
        fastq_file = r2_fastq[index],
        sample_id = sample_id,
        fastq_suffix = fastq_suffix
    }

    # if the index is passed, attach it to the bam file
    if (defined(i1_fastq)) {
      Array[File] non_optional_i1_fastq = select_first([i1_fastq])
      call Attach.Attach10xBarcodes as AttachBarcodes {
        input:
          r1_fastq = r1_fastq[index],
          i1_fastq = non_optional_i1_fastq[index],
          r2_unmapped_bam = FastqToUBam.bam_output,
          whitelist = whitelist,
          chemistry = chemistry
      }
    }

    # if the index is not passed, proceed without it.
    if (!defined(i1_fastq)) {
      call Attach.Attach10xBarcodes as AttachBarcodesNoIndex {
        input:
          r1_fastq = r1_fastq[index],
          r2_unmapped_bam = FastqToUBam.bam_output,
          whitelist = whitelist,
          chemistry = chemistry
      }
    }

    # This gets collected into an array outside of the scatter
    File barcoded_bam = select_first([AttachBarcodes.bam_output, AttachBarcodesNoIndex.bam_output])
  }

  scatter (bam in barcoded_bam) {
    call ScatterBam.ScatterBam as ScatterBamFiles {
      input:
        bam_to_scatter = bam,
        scatter_width = 32
    }

    Array[File] scattered_bams = ScatterBamFiles.scattered_bams
  }

  call ModifyGtf.ReplaceGeneNameWithGeneID as ModifyGtf {
    input:
      original_gtf = annotations_gtf
  }

  Array[File] flattened_scattered_bams = flatten(scattered_bams)

  call Split.SplitBamByCellBarcode {
    input:
      bams_to_split = flattened_scattered_bams
  }

  scatter (bam in SplitBamByCellBarcode.bam_output_array) {
    call StarAlignBam.StarAlignBamSingleEnd as StarAlign {
      input:
        bam_input = bam,
        tar_star_reference = tar_star_reference
    }

    if (counting_mode == "sc_rna") {
      call TagGeneExon.TagGeneExon as TagGenes {
        input:
          bam_input = StarAlign.bam_output,
          annotations_gtf = ModifyGtf.modified_gtf
      }
    }
    if (counting_mode == "sn_rna") {
      call TagGeneExon.TagReadWithGeneFunction as TagGeneFunction {
        input:
          bam_input = StarAlign.bam_output,
          annotations_gtf = ModifyGtf.modified_gtf,
          gene_name_tag = "GE",
          gene_strand_tag = "GS",
          gene_function_tag = "XF",
          use_strand_info = "false"
      }
    }

    call Picard.SortBamAndIndex as PreUMISort {
      input:
        bam_input = select_first([TagGenes.bam_output, TagGeneFunction.bam_output])
    }

    call UmiCorrection.CorrectUMItools as CorrectUMItools {
      input:
        bam_input = PreUMISort.bam_output,
        bam_index = PreUMISort.bam_index
    }

    call Picard.SortBamAndIndex as PreMergeSort {
      input:
        bam_input = CorrectUMItools.bam_output
    }

    call TagSortBam.GeneSortBam {
      input:
        bam_input = CorrectUMItools.bam_output
    }

    call TagSortBam.CellSortBam {
      input:
        bam_input = CorrectUMItools.bam_output
    }

    call Metrics.CalculateGeneMetrics {
      input:
        bam_input = GeneSortBam.bam_output
    }

    call Metrics.CalculateCellMetrics {
      input:
        bam_input = CellSortBam.bam_output,
        original_gtf = annotations_gtf
    }

    call Picard.SortBam as PreCountSort {
      input:
        bam_input = CorrectUMItools.bam_output,
        sort_order = "queryname"
    }

    call Count.CreateSparseCountMatrix {
      input:
        bam_input = PreCountSort.bam_output,
        gtf_file = ModifyGtf.modified_gtf
    }
  }

  call Merge.MergeSortBamFiles as MergeSorted {
    input:
      bam_inputs = PreMergeSort.bam_output,
      output_bam_filename = output_bam_basename + ".bam",
      sort_order = "coordinate"
  }

  call Metrics.MergeGeneMetrics {
    input:
      metric_files = CalculateGeneMetrics.gene_metrics
  }

  call Metrics.MergeCellMetrics {
    input:
      metric_files = CalculateCellMetrics.cell_metrics
  }

  call Count.MergeCountFiles {
    input:
      sparse_count_matrices = CreateSparseCountMatrix.sparse_count_matrix,
      row_indices = CreateSparseCountMatrix.row_index,
      col_indices = CreateSparseCountMatrix.col_index
  }

  call RunEmptyDrops.RunEmptyDrops {
    input:
      sparse_count_matrix = MergeCountFiles.sparse_count_matrix,
      row_index = MergeCountFiles.row_index,
      col_index = MergeCountFiles.col_index,
      emptydrops_lower = emptydrops_lower
  }

  call LoomUtils.OptimusLoomGeneration{
    input:
      sample_id = sample_id,
      annotation_file = annotations_gtf,
      cell_metrics = MergeCellMetrics.cell_metrics,
      gene_metrics = MergeGeneMetrics.gene_metrics,
      sparse_count_matrix = MergeCountFiles.sparse_count_matrix,
      cell_id = MergeCountFiles.row_index,
      gene_id = MergeCountFiles.col_index,
      empty_drops_result = RunEmptyDrops.empty_drops_result,
      counting_mode = counting_mode
  }

  output {
    # version of this pipeline
    String pipeline_version = version

    File bam = MergeSorted.output_bam
    File matrix = MergeCountFiles.sparse_count_matrix
    File matrix_row_index = MergeCountFiles.row_index
    File matrix_col_index = MergeCountFiles.col_index
    File cell_metrics = MergeCellMetrics.cell_metrics
    File gene_metrics = MergeGeneMetrics.gene_metrics
    File cell_calls = RunEmptyDrops.empty_drops_result

    # loom
    File loom_output_file = OptimusLoomGeneration.loom_output
  }
}
