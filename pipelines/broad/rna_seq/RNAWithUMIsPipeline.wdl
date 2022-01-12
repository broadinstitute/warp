version 1.0

import "../../../tasks/broad/UMIAwareDuplicateMarking.wdl" as UmiMD
import "../../../tasks/broad/RNAWithUMIsTasks.wdl" as tasks

## Copyright Broad Institute, 2021
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
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

workflow RNAWithUMIsPipeline {

	String pipeline_version = "0.1.0"

	input {
		File? bam
		File? r1_fastq
		File? r2_fastq
		String read1Structure
		String read2Structure
		File starIndex
		String output_basename
		File gtf

		# only needed if inputs are fastqs instead of ubam
		String? platform
		String? library_name
		String? platform_unit
		String? read_group_name
		String? sequencing_center = "BI"

		File ref
		File refIndex
		File refDict
		File refFlat
		File ribosomalIntervals
		File exonBedFile
	}

    call tasks.VerifyPipelineInputs {
        input:
            bam = bam,
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            library_name = library_name,
            platform = platform,
            platform_unit = platform_unit,
            read_group_name = read_group_name,
            sequencing_center = sequencing_center
    }

    if (VerifyPipelineInputs.fastq_run) {
        call tasks.FastqToUbam {
            input:
                r1_fastq = select_first([r1_fastq]),
                r2_fastq = select_first([r2_fastq]),
                bam_filename = output_basename,
                library_name = select_first([library_name]),
                platform = select_first([platform]),
                platform_unit = select_first([platform_unit]),
                read_group_name = select_first([read_group_name]),
                sequencing_center = select_first([sequencing_center])
        }
    }  

	File bam_to_use = select_first([bam, FastqToUbam.unmapped_bam])

	call tasks.ExtractUMIs {
		input:
			bam = bam_to_use,
			read1Structure = read1Structure,
			read2Structure = read2Structure
	}

	call tasks.STAR {
		input:
			bam = ExtractUMIs.bam_umis_extracted,
			starIndex = starIndex
	}

	call tasks.CopyReadGroupsToHeader {
		input:
			bam_with_readgroups = STAR.aligned_bam,
			bam_without_readgroups = STAR.transcriptome_bam
	}

	call UmiMD.UMIAwareDuplicateMarking {
		input:
			aligned_bam = STAR.aligned_bam,
			output_basename = output_basename
	}

	call UmiMD.UMIAwareDuplicateMarking as UMIAwareDuplicateMarkingTranscriptome {
		input:
			aligned_bam = CopyReadGroupsToHeader.output_bam,
			output_basename = output_basename + ".transcriptome"
	}

	### PLACEHOLDER for CROSSCHECK ###

	call tasks.GetSampleName {
		input:
			bam = bam_to_use
	}

	call tasks.rnaseqc2 {
		input:
			bam_file = UMIAwareDuplicateMarking.duplicate_marked_bam,
			genes_gtf = gtf,
			sample_id = GetSampleName.sample_name,
			exon_bed = exonBedFile
	}

	call tasks.CollectRNASeqMetrics {
		input:
			input_bam=UMIAwareDuplicateMarking.duplicate_marked_bam,
			input_bam_index=UMIAwareDuplicateMarking.duplicate_marked_bam_index,
			output_bam_prefix=GetSampleName.sample_name,
			ref_dict=refDict,
			ref_fasta=ref,
			ref_fasta_index=refIndex,
			ref_flat=refFlat,
			ribosomal_intervals=ribosomalIntervals,
	}

	call tasks.CollectMultipleMetrics {
		input:
			input_bam=UMIAwareDuplicateMarking.duplicate_marked_bam,
			input_bam_index=UMIAwareDuplicateMarking.duplicate_marked_bam_index,
			output_bam_prefix=GetSampleName.sample_name,
			ref_dict=refDict,
			ref_fasta=ref,
			ref_fasta_index=refIndex
	}

	# TODO: wire in fingerprint_summary_metrics once we have it. Using static example for now
	call tasks.MergeMetrics {
		input:
			alignment_summary_metrics=CollectMultipleMetrics.alignment_summary_metrics,
			insert_size_metrics=CollectMultipleMetrics.insert_size_metrics,
			picard_rna_metrics=CollectRNASeqMetrics.rna_metrics,
			duplicate_metrics=UMIAwareDuplicateMarking.duplicate_metrics,
			rnaseqc2_metrics=rnaseqc2.metrics,
			fingerprint_summary_metrics="gs://broad-gotc-test-storage/rna_seq/example.fingerprinting_summary_metrics",
			output_basename = GetSampleName.sample_name
	}


	output {
		File transcriptome_bam = UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam
		File transcriptome_bam_index = UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam_index
		File transcriptome_duplicate_metrics = UMIAwareDuplicateMarkingTranscriptome.duplicate_metrics
		File output_bam = UMIAwareDuplicateMarking.duplicate_marked_bam
		File output_bam_index = UMIAwareDuplicateMarking.duplicate_marked_bam_index
		File duplicate_metrics = UMIAwareDuplicateMarking.duplicate_metrics
		File rnaseqc2_gene_tpm = rnaseqc2.gene_tpm
		File rnaseqc2_gene_counts = rnaseqc2.gene_counts
		File rnaseqc2_exon_counts = rnaseqc2.exon_counts
		File rnaseqc2_fragment_size_histogram = rnaseqc2.fragment_size_histogram
		File rnaseqc2_metrics = rnaseqc2.metrics
		File picard_rna_metrics = CollectRNASeqMetrics.rna_metrics
		File picard_alignment_summary_metrics = CollectMultipleMetrics.alignment_summary_metrics
		File picard_insert_size_metrics = CollectMultipleMetrics.insert_size_metrics
		File picard_insert_size_histogram = CollectMultipleMetrics.insert_size_histogram
		File picard_base_distribution_by_cycle_metrics = CollectMultipleMetrics.base_distribution_by_cycle_metrics
		File picard_base_distribution_by_cycle_pdf = CollectMultipleMetrics.base_distribution_by_cycle_pdf
		File picard_quality_by_cycle_metrics = CollectMultipleMetrics.quality_by_cycle_metrics
		File picard_quality_by_cycle_pdf = CollectMultipleMetrics.quality_by_cycle_pdf
		File picard_quality_distribution_metrics = CollectMultipleMetrics.quality_distribution_metrics
		File picard_quality_distribution_pdf = CollectMultipleMetrics.quality_distribution_pdf
		File unified_metrics = MergeMetrics.unified_metrics
	}
}

