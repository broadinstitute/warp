version 1.0

import "../../../../../../tasks/broad/GermlineVariantDiscovery.wdl" as Calling
import "../../../../../../tasks/broad/Qc.wdl" as QC

workflow ReblockGVCF {

  String pipeline_version = "2.1.13"


  input {
    File gvcf
    File gvcf_index
    File? calling_interval_list
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? tree_score_cutoff
    String? annotations_to_keep_command
    String? annotations_to_remove_command
    Boolean? move_filters_to_genotypes
    String gvcf_file_extension = ".g.vcf.gz"
  }

  String gvcf_basename = basename(gvcf, gvcf_file_extension)

  call Calling.Reblock as Reblock {
    input:
      gvcf = gvcf,
      gvcf_index = gvcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      tree_score_cutoff = tree_score_cutoff,
      annotations_to_keep_command = annotations_to_keep_command,
      annotations_to_remove_command = annotations_to_remove_command,
      move_filters_to_genotypes = move_filters_to_genotypes,
      output_vcf_filename = gvcf_basename + ".rb.g.vcf.gz"
  }

    # Validate the (g)VCF output of HaplotypeCaller
    call QC.ValidateVCF as ValidateVCF {
      input:
        input_vcf = Reblock.output_vcf,
        input_vcf_index = Reblock.output_vcf_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        calling_interval_list = select_first([calling_interval_list, gvcf]), #nice trick so we don't have to pass around intervals; shouldn't be too much slower
        calling_interval_list_index = gvcf_index,
        calling_intervals_defined = defined(calling_interval_list),
        is_gvcf = true,
        extra_args = "--no-overlaps",
        gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

  output {
    File output_vcf = Reblock.output_vcf
    File output_vcf_index = Reblock.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}