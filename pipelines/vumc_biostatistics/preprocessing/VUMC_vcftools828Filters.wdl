version 1.0

#WORKFLOW DEFINITION
workflow VUMC_vcftools828Filters {
  input {
    File input_vcf
    String sample_name
    String vcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
  }
  
  call variantFilter {
    input: 
    input_vcf = input_vcf,
    sample_name = sample_name,
    vcftools_docker = vcftools_docker,
  }

  output {
    File output_828Filt_vcf = variantFilter.output_vcf
    File output_828Filt_vcf_index = variantFilter.output_vcf_index
  }
}

# TASK DEFINITIONS
# removing all indels and variants with no minor allele counts
task variantFilter{
  input{
    # Command parameters
    String sample_name
    File input_vcf
    
    # Runtime parameters
    String vcftools_docker
    Int machine_mem_gb = 4
    Int addtional_disk_space_gb = 50

  }

  Int disk_size = ceil(size(input_vcf, "GB")) + addtional_disk_space_gb
  String FilteredOutputVcf = "${sample_name}_828VariantFiltered.vcf.gz"
  String FilteredOutputVcfIndex = "${sample_name}_828VariantFiltered.vcf.gz.tbi"
  
  command <<<
    # remove the less than mac 1 and indels
    vcftools --gzvcf ~{input_vcf} --mac 1 --remove-indels --recode --stdout | bgzip -c > ~{FilteredOutputVcf}
    tabix ~{FilteredOutputVcf}
  >>>

  runtime{
    docker: vcftools_docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output{
    File output_vcf = "~{FilteredOutputVcf}"
    File output_vcf_index = "~{FilteredOutputVcfIndex}"
  }
}

