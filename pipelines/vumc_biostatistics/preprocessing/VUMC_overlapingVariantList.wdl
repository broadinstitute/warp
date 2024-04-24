version 1.0

#WORKFLOW DEFINITION
workflow VUMC_variantOverlaps {
  input {
    File first_input_vcf
    File second_input_vcf
    String first_sample_name
    String second_sample_name
    String vcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
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
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50

  }

  Int disk_size = 50
  String output_list = "${first_sample_name}_v_${second_sample_name}_shared_variants.txt"
  
  command <<<
    # compare and print out intersected variants
    tabix ~{first_input_vcf}
    tabix ~{second_input_vcf}
    bcftools isec -n +2 ~{first_input_vcf} ~{second_input_vcf} > ~{output_list}
  >>>

  runtime{
    docker: vcftools_docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output{
    File the_output_list = "~{output_list}"
  }
}