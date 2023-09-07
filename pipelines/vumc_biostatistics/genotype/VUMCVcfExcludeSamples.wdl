version 1.0

workflow VUMCVcfExcludeSamples {
  input {
    File input_vcf
    File input_vcf_index
    File exclude_samples
    String target_prefix
    String docker = "staphb/bcftools"

    String? project_id
    String? target_bucket
  }

  call BcftoolsExcludeSamples {
    input:
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      exclude_samples = exclude_samples,
      target_prefix = target_prefix,
      docker = docker
  }

  if(defined(target_bucket)){
    call CopyVcfFile {
      input:
        input_vcf = BcftoolsExcludeSamples.output_vcf,
        input_vcf_index = BcftoolsExcludeSamples.output_vcf_index,
        project_id = project_id,
        target_bucket = select_first([target_bucket])
    }
  }

  output {
    File output_vcf = select_first([CopyVcfFile.output_vcf, BcftoolsExcludeSamples.output_vcf])
    File output_vcf_index = select_first([CopyVcfFile.output_vcf_index, BcftoolsExcludeSamples.output_vcf_index])
  }
}

task BcftoolsExcludeSamples {
  input {
    File input_vcf
    File input_vcf_index
    File exclude_samples
    String target_prefix
    String docker = "staphb/bcftools"
  }

  Int disk_size = ceil(size(input_vcf, "GB") * 2) + 2

  command <<<

bcftools query -l ~{input_vcf} > all.id.txt

sort all.id.txt ~{exclude_samples} | uniq -d > remove.id.txt

bcftools view -S ^remove.id.txt -o ~{target_prefix}.vcf.gz ~{input_vcf}

bcftools index -t ~{target_prefix}.vcf.gz

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    File output_vcf = "~{target_prefix}.vcf.gz"
    File output_vcf_index = "~{target_prefix}.vcf.gz.tbi"
  }
}

task CopyVcfFile {
  input {
    String input_vcf
    String input_vcf_index

    String? project_id
    String target_bucket
  }

  String new_vcf = "~{target_bucket}/~{basename(input_vcf)}"
  String new_vcf_index = "~{target_bucket}/~{basename(input_vcf_index)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} cp ~{input_vcf} \
  ~{input_vcf_index} \
  ~{target_bucket}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String output_vcf = new_vcf
    String output_vcf_index = new_vcf_index
  }
}
