version 1.0

import "../../../tasks/skylab/StarAlign.wdl" as StarAlign
import "../../../tasks/skylab/FastqProcessing.wdl" as FastqProcessing
import "../../../tasks/skylab/Metrics.wdl" as Metrics
import "../../../tasks/skylab/LoomUtils.wdl" as LoomUtils
import "../../../tasks/skylab/CheckInputs.wdl" as OptimusInputChecks
import "../../../tasks/skylab/MergeSortBam.wdl" as Merge

## Copyright Broad Institute, 2022
##
## This WDL pipeline implements data processing for RNA with UMIs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. 


workflow SlideSeq {

    String pipeline_version = "2.1.3"

    input {
        Array[File] r1_fastq
        Array[File] r2_fastq
        Array[File]? i1_fastq
        String input_id
        String read_structure

        File tar_star_reference
        File annotations_gtf

        String output_bam_basename
        Boolean count_exons = true
        File bead_locations

    }

    parameter_meta {
        r1_fastq: "Array of Read 1 FASTQ files; forward read; contains cell barcodes and molecule barcodes"
        r2_fastq: "Array of Read 2 FASTQ files; reverse read; contains cDNA fragment generated from captured mRNA"
        i1_fastq: "Optional array of i1 FASTQ files; index read used for demultiplexing of multiple samples on one flow cell"
        input_id: "Name of sample matching this file; inserted into read group header"
        read_structure: "String used to specify the UMI (M) and Barcode (C) positions in the Read 1 FASTQ"
    }

    call StarAlign.STARGenomeRefVersion as ReferenceCheck {
    input:
      tar_star_reference = tar_star_reference
    }

    call Metrics.FastqMetricsSlideSeq as FastqMetrics {
        input:
            r1_fastq = r1_fastq,
            read_structure = read_structure,
            sample_id = input_id,
            whitelist = bead_locations
    }
    call FastqProcessing.FastqProcessingSlidSeq as SplitFastq {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            i1_fastq = i1_fastq,
            read_structure = read_structure,
            sample_id = input_id,
            whitelist = bead_locations
    }
    scatter(idx in range(length(SplitFastq.fastq_R1_output_array))) {
        call StarAlign.STARsoloFastqSlideSeq as STARsoloFastqSlideSeq {
            input:
                r1_fastq = [SplitFastq.fastq_R1_output_array[idx]],
                r2_fastq = [SplitFastq.fastq_R2_output_array[idx]],
                whitelist = bead_locations,
                tar_star_reference = tar_star_reference,
                output_bam_basename = output_bam_basename + "_" + idx,
                read_structure = read_structure,
                count_exons = count_exons
        }
    }
    call Merge.MergeSortBamFiles as MergeBam {
        input:
            bam_inputs = STARsoloFastqSlideSeq.bam_output,
            output_bam_filename = output_bam_basename + ".bam",
            sort_order = "coordinate"
    }
    call Metrics.CalculateGeneMetrics as GeneMetrics {
        input:
            bam_input = MergeBam.output_bam,
            input_id = input_id
    }
    call Metrics.CalculateUMIsMetrics as UMIsMetrics {
        input:
            bam_input = MergeBam.output_bam,
            input_id = input_id
    }

    call Metrics.CalculateCellMetrics as CellMetrics {
        input:
            bam_input = MergeBam.output_bam,
            original_gtf = annotations_gtf,
            input_id = input_id
    }

    call StarAlign.MergeStarOutput as MergeStarOutputs {
        input:
            barcodes = STARsoloFastqSlideSeq.barcodes,
            features = STARsoloFastqSlideSeq.features,
            matrix = STARsoloFastqSlideSeq.matrix,
            input_id = input_id
    }
    if ( !count_exons ) {
        call LoomUtils.OptimusLoomGeneration as SlideseqLoomGeneration{
            input:
                input_id = input_id,
                annotation_file = annotations_gtf,
                cell_metrics = CellMetrics.cell_metrics,
                gene_metrics = GeneMetrics.gene_metrics,
                sparse_count_matrix = MergeStarOutputs.sparse_counts,
                cell_id = MergeStarOutputs.row_index,
                gene_id = MergeStarOutputs.col_index,
                add_emptydrops_data = "no",
                pipeline_version = "SlideSeq_v~{pipeline_version}"
        }
    }
    if (count_exons) {
        call StarAlign.MergeStarOutput as MergeStarOutputsExons {
            input:
                barcodes = STARsoloFastqSlideSeq.barcodes_sn_rna,
                features = STARsoloFastqSlideSeq.features_sn_rna,
                matrix = STARsoloFastqSlideSeq.matrix_sn_rna,
                input_id = input_id
        }
        call LoomUtils.SingleNucleusOptimusLoomOutput as SlideseqLoomGenerationWithExons{
            input:
                input_id = input_id,
                annotation_file = annotations_gtf,
                cell_metrics = CellMetrics.cell_metrics,
                gene_metrics = GeneMetrics.gene_metrics,
                sparse_count_matrix = MergeStarOutputs.sparse_counts,
                cell_id = MergeStarOutputs.row_index,
                gene_id = MergeStarOutputs.col_index,
                sparse_count_matrix_exon = MergeStarOutputsExons.sparse_counts,
                cell_id_exon = MergeStarOutputsExons.row_index,
                gene_id_exon = MergeStarOutputsExons.col_index,
                pipeline_version = "SlideSeq_v~{pipeline_version}"
        }

    }

    File final_loom_output = select_first([SlideseqLoomGenerationWithExons.loom_output, SlideseqLoomGeneration.loom_output])

    output {
        String pipeline_version_out = pipeline_version
        File genomic_reference_version = ReferenceCheck.genomic_ref_version
        File bam = MergeBam.output_bam
        # sparse count matrix
        File matrix = MergeStarOutputs.sparse_counts
        File matrix_row_index = MergeStarOutputs.row_index
        File matrix_col_index = MergeStarOutputs.col_index

        File cell_metrics = CellMetrics.cell_metrics
        File gene_metrics = GeneMetrics.gene_metrics
        File umi_metrics =  UMIsMetrics.umi_metrics

        File fastq_barcode_distribution = FastqMetrics.barcode_distribution
        File fastq_umi_distribution = FastqMetrics.umi_distribution
        File fastq_reads_per_cell = FastqMetrics.numReads_perCell
        File fastq_reads_per_umi = FastqMetrics.numReads_perUMI


        # loom
        File? loom_output_file = final_loom_output

    }
}
