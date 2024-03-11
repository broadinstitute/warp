version 1.0

workflow VUMCVcfNumVariants {
  input {
    File input_vcf_index
    String docker = "staphb/bcftools"
  }

  call BcftoolsVcfNumVariants {
    input:
      input_vcf_index = input_vcf_index,
      docker = docker
  }

  output {
    Int num_variants = BcftoolsVcfNumVariants.num_variants
  }
}

task BcftoolsVcfNumVariants {
  input {
    File input_vcf_index
    String docker = "staphb/bcftools"
  }

  Int disk_size = ceil(size(input_vcf_index, "GB")) + 2
  String input_vcf_tmp = sub(input_vcf_index, ".csi$", "")
  String input_vcf = sub(input_vcf_tmp, ".tbi$", "")

  command <<<

vcf_index=$(echo ~{input_vcf_index} | sed -e 's/.csi$//g' | sed -e 's/.tbi$//g')

bcftools index -n $vcf_index > num_variants.txt

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    Int num_variants = read_int("num_variants.txt")
  }
}
