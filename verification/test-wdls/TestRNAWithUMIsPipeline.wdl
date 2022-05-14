version 1.0

import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../verification/VerifyRNAWithUMIs.wdl" as VerifyOptimus
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy
import "../../pipelines/broad/rna_seq/RNAWithUMIsPipeline.wdl" as RNAWithUMIsPipeline

workflow TestRNAWithUMIsPipeline {
  
  input {
      File? bam
      File? r1_fastq
      File? r2_fastq
      String read1Structure
      String read2Structure
      String output_basename
  
      String platform
      String library_name
      String platform_unit
      String read_group_name
      String sequencing_center = "BI"
  
      File starIndex
      File gtf
  
      File ref
      File refIndex
      File refDict
      File refFlat
      File ribosomalIntervals
      File exonBedFile
  
      File population_vcf
      File population_vcf_index
  
      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String vault_token_path
      String google_account_vault_path
  }

  meta {
    allowNestedInputs: true
  }

  call RNAWithUMIsPipeline.RNAWithUMIsPipeline {
    input:
        bam                  = bam,
        r1_fastq             = r1_fastq,
        r2_fastq             = r2_fastq,
        read1Structure       = read1Structure,
        read2Structure       = read2Structure,
        output_basename      = output_basename,
        platform             = platform,
        library_name         = library_name,
        platform_unit        = platform_unit,
        read_group_name      = read_group_name,
        sequencing_center    = sequencing_center,
        starIndex            = starIndex,
        gtf                  = gtf,
        ref                  = ref,
        refIndex             = refIndex,
        refDict              = refDict,
        refFlat              = refFlat,
        ribosomalIntervals   = ribosomalIntervals,
        exonBedFile          = exonBedFile,
        population_vcf       = population_vcf,
        population_vcf_index = population_vcf_index,
  }

  output{}
}