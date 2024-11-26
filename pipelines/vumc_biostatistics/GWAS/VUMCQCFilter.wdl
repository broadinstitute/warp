version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "../../../tasks/vumc_biostatistics/BioUtils.wdl" as BioUtils

workflow VUMCQCFilter {
  input {
    File input_pgen
    File input_pvar
    File input_psam

    String output_prefix

    String qc_option="--mac 100 --geno 0.1 --maf 0.1 --max-maf 0.9 --hwe 1e-15 --snps-only --not-chr 23-27"

    String? project_id
    String? target_gcp_folder
  }

  call BioUtils.PgenQCFilter {
    input:
      input_pgen = input_pgen,
      input_pvar = input_pvar,
      input_psam = input_psam,
      output_prefix = output_prefix,
      qc_option = qc_option
  }

  if(defined(target_gcp_folder)){
    String filtered_pgen = "~{PgenQCFilter.output_pgen}"
    String filtered_pvar = "~{PgenQCFilter.output_pvar}"
    String filtered_psam = "~{PgenQCFilter.output_psam}"

    call GcpUtils.MoveOrCopyThreeFiles as CopyFile {
      input:
        source_file1 = filtered_pgen,
        source_file2 = filtered_pvar,
        source_file3 = filtered_psam,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_pgen = select_first([CopyFile.output_file1, PgenQCFilter.output_pgen])
    File output_pvar = select_first([CopyFile.output_file2, PgenQCFilter.output_pvar])
    File output_psam = select_first([CopyFile.output_file3, PgenQCFilter.output_psam])
  }
}
