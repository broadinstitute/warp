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

  String pipeline_version = "1.0.2"

  input {
    File? bam
    File? r1_fastq
    File? r2_fastq
    String read1Structure
    String read2Structure
    String output_basename

    # The following inputs are only required if fastqs are given as input.
    String? platform
    String? library_name
    String? platform_unit
    String? read_group_name
    String? sequencing_center = "BI"

    File starIndex
    File gtf

    File ref
    File refIndex
    File refDict
    File refFlat
    File ribosomalIntervals
    File exonBedFile
  }

  parameter_meta {
    bam: "Read group-specific unmapped BAM file;  alternatively, paired-end FASTQ files (the `r1_fastq` and `r2_fastq` inputs) may be used"
    r1_fastq: "Read 1 FASTQ file; alternatively, the unmapped bam file (`bam` input) may be used as input"
    r2_fastq: "Read 2 FASTQ file; alternatively, the unmapped bam file (`bam` input) may be used as input"
    read1Structure: "String describing how the bases in a sequencing run should be allocated into logical reads for read 1"
    read2Structure: "String describing how the bases in a sequencing run should be allocated into logical reads for read 2"
    starIndex: "TAR file containing genome indices used for the STAR aligner"
    output_basename: "String used as a prefix in workflow output files"
    gtf: "Gene annotation file (GTF) used for the rnaseqc tool"
    platform: "String used to describe the sequencing platform; only required when using FASTQ files as input"
    library_name: "String used to describe the library; only required when using FASTQ files as input"
    platform_unit: "String used to describe the platform unit; only required when using FASTQ files as input"
    read_group_name: "String used to describe the read group name; only required when using FASTQ files as input"
    sequencing_center: "String used to describe the sequencing center; only required when using FASTQ files as input; default is set to 'BI'"
    ref: "FASTA file used for metric collection with Picard tools"
    refIndex: "FASTA index file used for metric collection with Picard tools"
    refDict: "Dictionary file used for metric collection with Picard tools"
    refFlat: "refFlat file used for metric collection with Picard tools"
    ribosomalIntervals: "Intervals file used for RNA metric collection with Picard tools"
    exonBedFile: "Bed file used for fragment size calculations in the rnaseqc tool; contains non-overlapping exons"
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
      input_bam = UMIAwareDuplicateMarking.duplicate_marked_bam,
      input_bam_index = UMIAwareDuplicateMarking.duplicate_marked_bam_index,
      output_bam_prefix = GetSampleName.sample_name,
      ref_dict = refDict,
      ref_fasta = ref,
      ref_fasta_index = refIndex,
      ref_flat = refFlat,
      ribosomal_intervals = ribosomalIntervals,
  }

  call tasks.CollectMultipleMetrics {
    input:
      input_bam = UMIAwareDuplicateMarking.duplicate_marked_bam,
      input_bam_index = UMIAwareDuplicateMarking.duplicate_marked_bam_index,
      output_bam_prefix = GetSampleName.sample_name,
      ref_dict = refDict,
      ref_fasta = ref,
      ref_fasta_index = refIndex
  }

  output {
    String sample_name = GetSampleName.sample_name
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
  }
}

