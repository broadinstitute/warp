version 1.0

import "./Utils.wdl" as Utils

workflow VUMCVcfExcludeSamples {
  input {
    File input_vcf
    File exclude_samples
    String target_prefix
    String target_suffix = ".vcf.gz"
    String docker = "staphb/bcftools"

    String? project_id
    String? target_gcp_folder
  }

  call BcftoolsExcludeSamples {
    input:
      input_vcf = input_vcf,
      exclude_samples = exclude_samples,
      target_prefix = target_prefix,
      target_suffix = target_suffix,
      docker = docker
  }

  if(defined(target_gcp_folder)){
    call Utils.MoveOrCopyTwoFiles as CopyVCF {
      input:
        source_file1 = BcftoolsExcludeSamples.output_vcf,
        source_file2 = BcftoolsExcludeSamples.output_vcf_index,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder]),
    }
  }

  output {
    File output_vcf = select_first([CopyVCF.output_file1, BcftoolsExcludeSamples.output_vcf])
    File output_vcf_index = select_first([CopyVCF.output_file2, BcftoolsExcludeSamples.output_vcf_index])
  }
}

task BcftoolsExcludeSamples {
  input {
    File input_vcf
    File exclude_samples
    String target_prefix
    String target_suffix = ".vcf.bgz"
    String docker = "staphb/bcftools"
  }

  Int disk_size = ceil(size(input_vcf, "GB") * 2) + 2
  String new_vcf = target_prefix + target_suffix

  command <<<

bcftools query -l ~{input_vcf} > all.id.txt

grep -Fxf all.id.txt ~{exclude_samples} | sort | uniq > remove.id.txt

bcftools view -S ^remove.id.txt -o ~{new_vcf} ~{input_vcf}

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

