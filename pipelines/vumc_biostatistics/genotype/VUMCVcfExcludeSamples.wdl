version 1.0

workflow VUMCVcfExcludeSamples {
  input {
    File input_vcf
    File input_vcf_index
    File exclude_samples
    String target_prefix
    String docker = "staphb/bcftools"
  }

  call BcftoolsExcludeSamples {
    input:
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      exclude_samples = exclude_samples,
      target_prefix = target_prefix,
      docker = docker
  }

  output {
    File output_vcf = BcftoolsExcludeSamples.output_vcf
    File output_vcf_index = BcftoolsExcludeSamples.output_vcf_index
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