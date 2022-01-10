version 1.0

import "RNAWithUMIsPipeline.wdl" as RNAWithUMIsPipeline
import "../../../tasks/broad/Utilities.wdl" as utils


workflow BroadInternalRNAWithUMIsPipeline {

  String pipeline_version = "0.1.0"

  input {
    # input needs to be either "hg19" or "b38"
    String reference_build

    # RNAWithUMIsPipeline inputs
    File? bam
    File? r1_fastq
    File? r2_fastq
    String read1Structure
    String read2Structure
    String output_basename

    # only needed if inputs are fastqs instead of ubam
    String? platform
    String? library_name
    String? platform_unit
    String? read_group_name
    String? sequencing_center = "BI"

  }

  parameter_meta {
    bam: "Read group-specific unmapped BAM file;  alternatively, paired-end FASTQ files (the `r1_fastq` and `r2_fastq` inputs) may be used"
    r1_fastq: "Read 1 FASTQ file; alternatively, the unmapped bam file (`bam` input) may be used as input"
    r2_fastq: "Read 2 FASTQ file; alternatively, the unmapped bam file (`bam` input) may be used as input"
    read1Structure: "String describing how the bases in a sequencing run should be allocated into logical reads for read 1"
    read2Structure: "String describing how the bases in a sequencing run should be allocated into logical reads for read 2"
    output_basename: "String used as a prefix in workflow output files"
    platform: "String used to describe the sequencing platform; only required when using FASTQ files as input"
    library_name: "String used to describe the library; only required when using FASTQ files as input"
    platform_unit: "String used to describe the platform unit; only required when using FASTQ files as input"
    read_group_name: "String used to describe the read group name; only required when using FASTQ files as input"
    sequencing_center: "String used to describe the sequencing center; only required when using FASTQ files as input; default is set to “BI”"
  }

  # make sure either hg19 or b38 is supplied as reference_build input
  if (reference_build != "hg19") {
    if (reference_build != "b38") {
      call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
        input:
          message = "reference_build must be supplied with either 'hg19' or 'b38'."
      }
    }
  }

 # if reference_build=hg19, use hg19 references:
 if (reference_build == "hg19") {
   call RNAWithUMIsPipeline.RNAWithUMIsPipeline as RNAWithUMIsPipelineHg19 {
     input:
       bam = bam,
       r1_fastq = r1_fastq,
       r2_fastq = r2_fastq,
       read1Structure = read1Structure,
       read2Structure = read2Structure,
       starIndex = "gs://broad-gotc-test-storage/rna_seq/hg19/STAR_genome_hg19_v19.tar.gz",
       output_basename = output_basename,
       gtf = "gs://broad-gotc-test-storage/rna_seq/hg19/gencode.v19.genes.v7.collapsed_only.patched_contigs.gtf",
       platform = platform,
       library_name = library_name,
       platform_unit = platform_unit,
       read_group_name = read_group_name,
       sequencing_center = "BI",
       ref = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta",
       refIndex = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
       refDict = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict",
       refFlat = "gs://broad-gotc-test-storage/rna_seq/hg19/Homo_sapiens_assembly19.refFlat.txt",
       ribosomalIntervals = "gs://broad-gotc-test-storage/rna_seq/hg19/Homo_sapiens_assembly19.rRNA.interval_list",
       exonBedFile = "gs://broad-gotc-test-storage/rna_seq/hg19/gencode.v19.hg19.insert_size_intervals_geq1000bp.bed"
   }
 }

 # if reference_build=b38, use hg38 refernces:
 if (reference_build == "b38") {
   call RNAWithUMIsPipeline.RNAWithUMIsPipeline as RNAWithUMIsPipelineHg38 {
     input:
       bam = bam,
       r1_fastq = r1_fastq,
       r2_fastq = r2_fastq,
       read1Structure = read1Structure,
       read2Structure = read2Structure,
       starIndex = "gs://broad-gotc-test-storage/rna_seq/hg38/STAR_genome_GRCh38_noALT_noHLA_noDecoy_v26_oh149.tar.gz",
       output_basename = output_basename,
       gtf = "gs://broad-gotc-test-storage/rna_seq/hg38/gencode.v26.GRCh38.genes.collapsed_only.gtf",
       platform = platform,
       library_name = library_name,
       platform_unit = platform_unit,
       read_group_name = read_group_name,
       sequencing_center = "BI",
       ref = "gs://broad-gotc-test-storage/rna_seq/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta",
       refIndex = "gs://broad-gotc-test-storage/rna_seq/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta.fai",
       refDict = "gs://broad-gotc-test-storage/rna_seq/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.dict",
       refFlat = "gs://broad-gotc-test-storage/rna_seq/hg38/GRCh38_gencode.v27.refFlat.txt",
       ribosomalIntervals = "gs://broad-gotc-test-storage/rna_seq/hg38/gencode.v26.rRNA.withMT.interval_list",
       exonBedFile = "gs://broad-gotc-test-storage/rna_seq/hg38/gencode.v26.GRCh38.insert_size_intervals_geq1000bp.bed"
   }
 }
}
