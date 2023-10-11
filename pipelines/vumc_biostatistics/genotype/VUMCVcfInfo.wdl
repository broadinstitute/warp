version 1.0

workflow VUMCVcfInfo {
  input {
    File input_vcf
    File input_vcf_index
    String docker = "staphb/bcftools"
  }

  call BcftoolsVcfInfo {
    input:
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      docker = docker
  }

  output {
    File sample_file = BcftoolsVcfInfo.sample_file
    Int num_samples = BcftoolsVcfInfo.num_samples
    Int num_variants = BcftoolsVcfInfo.num_variants
  }
}

task BcftoolsVcfInfo {
  input {
    File input_vcf
    File input_vcf_index
    String docker = "staphb/bcftools"
  }

  Int disk_size = ceil(size(input_vcf, "GB")) + 2
  String output_sample_file = basename(input_vcf) + ".samples.txt"

  command <<<

bcftools query -l ~{input_vcf} > ~{output_sample_file}

cat ~{output_sample_file} | wc -l > num_samples.txt

bcftools index -n ~{input_vcf} > num_variants.txt

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    File sample_file = "~{output_sample_file}"
    Int num_samples = read_int("num_samples.txt")
    Int num_variants = read_int("num_variants.txt")
  }
}
