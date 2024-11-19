version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "../agd/AgdUtils.wdl" as AgdUtils

workflow VUMCAnnovar {
  input {
    File input_vcf

    String target_prefix

    String? billing_project_id
    String? target_gcp_folder
  }

  call AgdUtils.Annovar {
    input:
      input_vcf = input_vcf,
      target_prefix = target_prefix
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyOneFile as CopyFile {
      input:
        source_file = Annovar.annovar_file,
        is_move_file = false,
        project_id = billing_project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File annovar_file = select_first([CopyFile.output_file, Annovar.annovar_file])
  }
}
