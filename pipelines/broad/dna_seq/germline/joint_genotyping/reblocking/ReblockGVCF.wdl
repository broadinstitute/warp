version 1.0

import "../../../../../../tasks/broad/GermlineVariantDiscovery.wdl" as Calling
import "../../../../../../tasks/broad/QC.wdl" as QC

workflow ReblockGVCF {

  String pipeline_version = "1.1.1"

  input {
    File gvcf
    File gvcf_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
  }

  String gvcf_basename = basename(gvcf, ".g.vcf.gz")

  call Calling.Reblock as Reblock {
    input:
      gvcf = gvcf,
      gvcf_index = gvcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      output_vcf_filename = gvcf_basename + ".reblocked.g.vcf.gz"
  }

    # Validate the (g)VCF output of HaplotypeCaller
    call QC.ValidateVCF as ValidateVCF {
      input:
        input_vcf = Reblock.output_vcf,
        input_vcf_index = Reblock.output_vcf_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        calling_interval_list = gvcf, #nice trick so we don't have to pass around intervals; shouldn't be too much slower
        is_gvcf = true,
        extra_args = "--no-overlaps",
        gatk_docker = "us.gcr.io/broad-dsde-methods/validate_reblocking@sha256:ede2cae0d345c6ad2690db99d1d6485adc2b52f98deb5e821bd64baf3df63602"
    }

  output {
    File output_vcf = Reblock.output_vcf
    File output_vcf_index = Reblock.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
