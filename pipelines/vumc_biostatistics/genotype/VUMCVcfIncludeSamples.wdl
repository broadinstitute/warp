version 1.0

import "./Utils.wdl" as Utils

workflow VUMCVcfIncludeSamples {
  input {
    File input_vcf
    File input_vcf_index
    File include_samples
    String target_prefix
    String target_suffix = ".vcf.gz"
    String docker = "staphb/bcftools"

    String? project_id
    String? target_bucket
    String? genoset
  }

  call BcftoolsIncludeSamples {
    input:
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      include_samples = include_samples,
      target_prefix = target_prefix,
      target_suffix = target_suffix,
      docker = docker
  }

  if(defined(target_bucket)){
    call Utils.CopyVcfFile {
      input:
        input_vcf = BcftoolsIncludeSamples.output_vcf,
        input_vcf_index = BcftoolsIncludeSamples.output_vcf_index,
        project_id = project_id,
        target_bucket = select_first([target_bucket]),
        genoset = select_first([genoset]),
    }
  }

  output {
    File output_vcf = select_first([CopyVcfFile.output_vcf, BcftoolsIncludeSamples.output_vcf])
    File output_vcf_index = select_first([CopyVcfFile.output_vcf_index, BcftoolsIncludeSamples.output_vcf_index])
  }
}

task BcftoolsIncludeSamples {
  input {
    File input_vcf
    File input_vcf_index
    File include_samples
    String target_prefix
    String target_suffix = ".vcf.gz"
    String docker = "staphb/bcftools"
  }

  Int disk_size = ceil(size(input_vcf, "GB") * 2) + 2
  String new_vcf = target_prefix + target_suffix

  command <<<

bcftools query -l ~{input_vcf} > all.id.txt

sort all.id.txt ~{include_samples} | uniq -d > keep.id.txt

bcftools view -S keep.id.txt -o ~{new_vcf} ~{input_vcf}

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
