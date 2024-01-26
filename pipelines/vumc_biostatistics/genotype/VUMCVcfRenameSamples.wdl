version 1.0

import "./Utils.wdl" as Utils

workflow VUMCVcfRenameSamples {
  input {
    File input_vcf
    File replace_samples
    String target_prefix
    String target_suffix = ".vcf.gz"
    String docker = "staphb/bcftools"

    String? project_id
    String? target_bucket
    String? genoset
  }

  call BcftoolsReplaceHeader {
    input:
      input_vcf = input_vcf,
      replace_samples = replace_samples,
      target_prefix = target_prefix,
      target_suffix = target_suffix,
      docker = docker
  }

  if(defined(target_bucket)){
    call Utils.MoveOrCopyVcfFile {
      input:
        input_vcf = BcftoolsReplaceHeader.output_vcf,
        input_vcf_index = BcftoolsReplaceHeader.output_vcf_index,
        is_move_file = true,
        project_id = project_id,
        target_bucket = select_first([target_bucket]),
        genoset = select_first([genoset]),
    }
  }

  output {
    File output_vcf = select_first([MoveOrCopyVcfFile.output_vcf, BcftoolsReplaceHeader.output_vcf])
    File output_vcf_index = select_first([MoveOrCopyVcfFile.output_vcf_index, BcftoolsReplaceHeader.output_vcf_index])
  }
}

task BcftoolsReplaceHeader {
  input {
    File input_vcf
    File replace_samples
    String target_prefix
    String target_suffix = ".vcf.gz"
    String docker = "staphb/bcftools"
    Float disk_factor = 3.0
  }

  Int disk_size = ceil(size(input_vcf, "GB") * disk_factor) + 2
  String new_vcf = target_prefix + target_suffix

  command <<<

echo bcftools reheader -s ~{replace_samples} -o ~{new_vcf} ~{input_vcf}
bcftools reheader -s ~{replace_samples} -o ~{new_vcf} ~{input_vcf}

echo "build index"
bcftools index -t ~{new_vcf}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    File output_vcf = "~{new_vcf}"
    File output_vcf_index = "~{new_vcf}.tbi"
  }
}
