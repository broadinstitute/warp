version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils
import "../../../tasks/vumc_biostatistics/BioUtils.wdl" as BioUtils

# QC filter pgen files and output samples and variants list
workflow VUMCQCFilterList {
  input {
    File input_pgen
    File input_pvar
    File input_psam

    String output_prefix

    String qc_option="--mac 100 --geno 0.1 --maf 0.1 --max-maf 0.9 --hwe 1e-15 --snps-only --not-chr 23-27"

    String? project_id
    String? target_gcp_folder
  }

  call BioUtils.PgenQCFilterList {
    input:
      input_pgen = input_pgen,
      input_pvar = input_pvar,
      input_psam = input_psam,
      output_prefix = output_prefix,
      qc_option = qc_option
  }

  if(defined(target_gcp_folder)){
    String snplist_file = "~{PgenQCFilterList.output_snplist}"
    String samples_file = "~{PgenQCFilterList.output_samples}"

    call GcpUtils.MoveOrCopyTwoFiles as CopyFile {
      input:
        source_file1 = snplist_file,
        source_file2 = samples_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File qc_filtered_snplist = select_first([CopyFile.output_file1, PgenQCFilterList.output_snplist])
    File qc_filtered_samples = select_first([CopyFile.output_file2, PgenQCFilterList.output_samples])
  }
}
