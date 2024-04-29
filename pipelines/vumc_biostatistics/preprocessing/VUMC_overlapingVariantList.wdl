version 1.0

#WORKFLOW DEFINITION
workflow VUMC_variantOverlaps {
  input {
    File first_input_vcf
    File second_input_vcf
    String first_sample_name
    String second_sample_name
    String vcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int addtional_disk_space_gb = 50
    Int machine_mem_gb = 4
  }
  
  call varOverlap {
    input: 
    first_input_vcf = first_input_vcf,
    second_input_vcf = second_input_vcf,
    first_sample_name = first_sample_name,
    second_sample_name = second_sample_name,
    vcftools_docker = vcftools_docker,
  }

  output {
    File shared_list = varOverlap.the_output_list
    File output_vcf1 = varOverlap.first_output_vcf
    File output_vcfIndex1 = varOverlap.first_output_vcfIndex
    File output_vcf2 = varOverlap.second_output_vcf
    File output_vcfIndex2 = varOverlap.second_output_vcfIndex
    
  }
}

# TASK DEFINITIONS
# removing all indels and variants with no minor allele counts
task varOverlap{
  input{
    # Command parameters
    String first_sample_name
    String second_sample_name
    File first_input_vcf
    File second_input_vcf
    
    # Runtime parameters
    String vcftools_docker
    Int machine_mem_gb
    Int addtional_disk_space_gb

  }

  Int disk_size = ceil(size(first_input_vcf, "GB")) + ceil(size(second_input_vcf, "GB")) + addtional_disk_space_gb
  String output_list = "${first_sample_name}_v_${second_sample_name}_shared_variants.txt"
  String output_first_vcf = "${first_sample_name}_shared_variants.vcf.gz"
  String output_first_vcfIndex = "${first_sample_name}_shared_variants.vcf.gz.tbi"
  String output_second_vcf = "${second_sample_name}_shared_variants.vcf.gz"
  String output_second_vcfIndex = "${second_sample_name}_shared_variants.vcf.gz.tbi"
  
  command <<<
    # compare and print out intersected variants
    tabix ~{first_input_vcf}
    tabix ~{second_input_vcf}
    bcftools isec -n =2 ~{first_input_vcf} ~{second_input_vcf} > ~{output_list}
    cut -f1,2 ~{output_list}>tmp.txt
    vcftools --gzvcf ~{first_input_vcf} --positions tmp.txt --recode --stdout |bgzip -c>~{output_first_vcf}
    vcftools --gzvcf ~{second_input_vcf} --positions tmp.txt --recode --stdout |bgzip -c>~{output_second_vcf}
    tabix ~{output_first_vcf}
    tabix ~{output_second_vcf}
  >>>

  runtime{
    docker: vcftools_docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output{
    File the_output_list = "~{output_list}"
    File first_output_vcf = "~{output_first_vcf}"
    File first_output_vcfIndex = "~{output_first_vcf}.tbi"
    File second_output_vcf = "~{output_second_vcf}"
    File second_output_vcfIndex = "~{output_second_vcf}.tbi"
  }
}