version 1.0

import "../../../pipelines/broad/rna_seq/RNAWithUMIsPipeline.wdl" as RNAWithUMIsPipeline
import "../../../pipelines/broad/qc/CheckFingerprint.wdl" as FP
import "../../../tasks/broad/RNAWithUMIsTasks.wdl" as tasks
import "../../../tasks/broad/Utilities.wdl" as utils

workflow BroadInternalRNAWithUMIsPipeline {

  String pipeline_version = "0.1.0"

  input {
    # input needs to be either "hg19" or "hg38"
    String reference_build

    String sample_alias
    String sample_lsid

    # RNAWithUMIsPipeline inputs
    File r1_fastq
    File r2_fastq
    String read1Structure
    String read2Structure
    String output_basename

    String platform
    String library_name
    String platform_unit
    String read_group_name
    String sequencing_center = "BI"

    String environment
    File vault_token_path
  }

  File starIndex = if (reference_build == "hg19") then "gs://broad-gotc-test-storage/rna_seq/hg19/STAR_genome_hg19_v19.tar.gz" else "gs://broad-gotc-test-storage/rna_seq/hg38/STAR_genome_GRCh38_noALT_noHLA_noDecoy_v26_oh149.tar.gz"
  File gtf = if (reference_build == "hg19") then "gs://broad-gotc-test-storage/rna_seq/hg19/gencode.v19.genes.v7.collapsed_only.patched_contigs.gtf" else "gs://broad-gotc-test-storage/rna_seq/hg38/gencode.v26.GRCh38.genes.collapsed_only.gtf"
  File ref = if (reference_build == "hg19") then "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta" else "gs://broad-gotc-test-storage/rna_seq/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"
  File refIndex = if (reference_build == "hg19") then "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai" else "gs://broad-gotc-test-storage/rna_seq/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai"
  File refDict = if (reference_build == "hg19") then "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict" else "gs://broad-gotc-test-storage/rna_seq/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict"
  File refFlat = if (reference_build == "hg19") then "gs://broad-gotc-test-storage/rna_seq/hg19/Homo_sapiens_assembly19.refFlat.txt" else "gs://broad-gotc-test-storage/rna_seq/hg38/GRCh38_gencode.v27.refFlat.txt"
  File haplotype_database_file = if (reference_build == "hg19") then "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.haplotype_database.txt" else "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt"
  File ribosomalIntervals = if (reference_build == "hg19") then "gs://broad-gotc-test-storage/rna_seq/hg19/Homo_sapiens_assembly19.rRNA.interval_list" else "gs://broad-gotc-test-storage/rna_seq/hg38/gencode.v26.rRNA.withMT.interval_list"
  File exonBedFile = if (reference_build == "hg19") then "gs://broad-gotc-test-storage/rna_seq/hg19/gencode.v19.hg19.insert_size_intervals_geq1000bp.bed" else "gs://broad-gotc-test-storage/rna_seq/hg38/gencode.v26.GRCh38.insert_size_intervals_geq1000bp.bed"

  parameter_meta {
    reference_build: "String used to define the reference genome build; should be set to 'hg19' or 'hg38'"
    sample_alias: "The sample alias (aka the Collaborator Sample ID) the name of the sample"
    sample_lsid: "The sample lsid (an identifier used to retrieve fingerrints from Mercury)"
    r1_fastq: "Read 1 FASTQ file"
    r2_fastq: "Read 2 FASTQ file"
    read1Structure: "String describing how the bases in a sequencing run should be allocated into logical reads for read 1"
    read2Structure: "String describing how the bases in a sequencing run should be allocated into logical reads for read 2"
    output_basename: "String used as a prefix in workflow output files"
    platform: "String used to describe the sequencing platform; only required when using FASTQ files as input"
    library_name: "String used to describe the library; only required when using FASTQ files as input"
    platform_unit: "String used to describe the platform unit; only required when using FASTQ files as input"
    read_group_name: "String used to describe the read group name; only required when using FASTQ files as input"
    sequencing_center: "String used to describe the sequencing center; only required when using FASTQ files as input; default is set to 'BI'"
    environment: "The environment (dev or prod) used for determining which service to use to retrieve Mercury fingerprints"
    vault_token_path: "The path to the vault token used for accessing the Mercury Fingerprint Store"
  }

  # make sure either hg19 or hg38 is supplied as reference_build input
  if ((reference_build != "hg19") && (reference_build != "hg38")) {
    call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
      input:
        message = "reference_build must be supplied with either 'hg19' or 'hg38'."
    }
  }

  call RNAWithUMIsPipeline.RNAWithUMIsPipeline as RNAWithUMIsPipeline {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      read1Structure = read1Structure,
      read2Structure = read2Structure,
      starIndex = starIndex,
      output_basename = output_basename,
      gtf = gtf,
      platform = platform,
      library_name = library_name,
      platform_unit = platform_unit,
      read_group_name = read_group_name,
      sequencing_center = sequencing_center,
      ref = ref,
      refIndex = refIndex,
      refDict = refDict,
      refFlat = refFlat,
      ribosomalIntervals = ribosomalIntervals,
      exonBedFile = exonBedFile
  }

  call FP.CheckFingerprint as CheckFingerprint {
    input:
      input_bam = RNAWithUMIsPipeline.output_bam,
      input_bam_index = RNAWithUMIsPipeline.output_bam_index,
      sample_alias = sample_alias,
      sample_lsid = sample_lsid,
      output_basename = output_basename,
      ref_fasta = ref,
      ref_fasta_index = refIndex,
      ref_dict = refDict,
      read_fingerprint_from_mercury = true,
      haplotype_database_file = haplotype_database_file,
      environment = environment,
      vault_token_path = vault_token_path
  }

  call tasks.MergeMetrics {
    input:
      alignment_summary_metrics = RNAWithUMIsPipeline.picard_alignment_summary_metrics,
      insert_size_metrics = RNAWithUMIsPipeline.picard_insert_size_metrics,
      picard_rna_metrics = RNAWithUMIsPipeline.picard_rna_metrics,
      duplicate_metrics = RNAWithUMIsPipeline.duplicate_metrics,
      rnaseqc2_metrics = RNAWithUMIsPipeline.rnaseqc2_metrics,
      fingerprint_summary_metrics = CheckFingerprint.fingerprint_summary_metrics_file,
      output_basename = RNAWithUMIsPipeline.sample_name
  }

  output {
    File transcriptome_bam = RNAWithUMIsPipeline.transcriptome_bam
    File transcriptome_bam_index = RNAWithUMIsPipeline.transcriptome_bam_index
    File transcriptome_duplicate_metrics = RNAWithUMIsPipeline.transcriptome_duplicate_metrics
    File output_bam = RNAWithUMIsPipeline.output_bam
    File output_bam_index = RNAWithUMIsPipeline.output_bam_index
    File duplicate_metrics = RNAWithUMIsPipeline.duplicate_metrics
    File rnaseqc2_gene_tpm = RNAWithUMIsPipeline.rnaseqc2_gene_tpm
    File rnaseqc2_gene_counts = RNAWithUMIsPipeline.rnaseqc2_gene_counts
    File rnaseqc2_exon_counts = RNAWithUMIsPipeline.rnaseqc2_exon_counts
    File rnaseqc2_fragment_size_histogram = RNAWithUMIsPipeline.rnaseqc2_fragment_size_histogram
    File rnaseqc2_metrics = RNAWithUMIsPipeline.rnaseqc2_metrics
    File picard_rna_metrics = RNAWithUMIsPipeline.picard_rna_metrics
    File picard_alignment_summary_metrics = RNAWithUMIsPipeline.picard_alignment_summary_metrics
    File picard_insert_size_metrics = RNAWithUMIsPipeline.picard_insert_size_metrics
    File picard_insert_size_histogram = RNAWithUMIsPipeline.picard_insert_size_histogram
    File picard_base_distribution_by_cycle_metrics = RNAWithUMIsPipeline.picard_base_distribution_by_cycle_metrics
    File picard_base_distribution_by_cycle_pdf = RNAWithUMIsPipeline.picard_base_distribution_by_cycle_pdf
    File picard_quality_by_cycle_metrics = RNAWithUMIsPipeline.picard_quality_by_cycle_metrics
    File picard_quality_by_cycle_pdf = RNAWithUMIsPipeline.picard_quality_by_cycle_pdf
    File picard_quality_distribution_metrics = RNAWithUMIsPipeline.picard_quality_distribution_metrics
    File picard_quality_distribution_pdf = RNAWithUMIsPipeline.picard_quality_distribution_pdf
    File unified_metrics = MergeMetrics.unified_metrics
  }
}
