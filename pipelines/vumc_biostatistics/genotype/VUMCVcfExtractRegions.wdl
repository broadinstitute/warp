version 1.0

workflow VUMCVcfExtractRegions {
  input {
    File input_vcf
    File input_vcf_index
    File include_region_bed

    String target_prefix
    String target_suffix = ".vcf.gz"

    String? project_id
    String? target_gcp_folder
  }

  if(!defined(target_gcp_folder)){
    call BcftoolsExtractRegions {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        include_region_bed = include_region_bed,
        target_prefix = target_prefix,
        target_suffix = target_suffix
    }
  }

  if(defined(target_gcp_folder)){
    call BcftoolsExtractRegionsGcp {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        include_region_bed = include_region_bed,
        target_prefix = target_prefix,
        target_suffix = target_suffix,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_vcf = select_first([BcftoolsExtractRegions.output_vcf, BcftoolsExtractRegionsGcp.output_vcf])
    File output_vcf_index = select_first([BcftoolsExtractRegions.output_vcf_index, BcftoolsExtractRegionsGcp.output_vcf_index])
    File output_vcf_sample = select_first([BcftoolsExtractRegions.output_vcf_sample, BcftoolsExtractRegionsGcp.output_vcf_sample])
    Int output_vcf_num_samples = select_first([BcftoolsExtractRegions.output_vcf_num_samples, BcftoolsExtractRegionsGcp.output_vcf_num_samples])
    Int output_vcf_num_variants = select_first([BcftoolsExtractRegions.output_vcf_num_variants, BcftoolsExtractRegionsGcp.output_vcf_num_variants])
  }
}

task BcftoolsExtractRegions {
  input {
    File input_vcf
    File input_vcf_index
    File include_region_bed

    String target_prefix
    String target_suffix
    
    String docker = "shengqh/hail_gcp:20240213"
    Float disk_factor = 3.0
    Int preemptible = 1
    Int cpu = 8
  }

  Int disk_size = ceil(size(input_vcf, "GB") * disk_factor) + 2
  Int memory_gb = 2 * cpu

  String target_vcf = target_prefix + target_suffix
  String target_vcf_index = target_vcf + ".tbi"
  String target_sample_file = target_vcf + ".samples.txt"

  command <<<

echo bcftools filter -R ~{include_region_bed} --threads ~{cpu} -o ~{target_vcf} ~{input_vcf}
bcftools filter -R ~{include_region_bed} --threads ~{cpu} -o ~{target_vcf} ~{input_vcf}

echo "build index"
bcftools index -t --threads ~{cpu} ~{target_vcf}

bcftools query -l ~{target_vcf} > ~{target_sample_file}

cat ~{target_sample_file} | wc -l > num_samples.txt

bcftools index -n ~{target_vcf} > num_variants.txt

>>>

  runtime {
    docker: docker
    preemptible: preemptible
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    # The output has to be defined as File, otherwise the file would not be delocalized
    File output_vcf = "~{target_vcf}"
    File output_vcf_index = "~{target_vcf_index}"
    File output_vcf_sample = "~{target_sample_file}"
    Int output_vcf_num_samples = read_int("num_samples.txt")
    Int output_vcf_num_variants = read_int("num_variants.txt")
  }
}

task BcftoolsExtractRegionsGcp {
  input {
    File input_vcf
    File input_vcf_index
    File include_region_bed

    String target_prefix
    String target_suffix
    
    String? project_id
    String target_gcp_folder

    String docker = "shengqh/hail_gcp:20240213"
    Float disk_factor = 3.0
    Int preemptible = 1
    Int cpu = 8
  }

  Int disk_size = ceil(size(input_vcf, "GB") * disk_factor) + 2
  Int memory_gb = 2 * cpu

  String target_vcf = target_prefix + target_suffix
  String target_vcf_index = target_vcf + ".tbi"
  String target_sample_file = target_vcf + ".samples.txt"

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")
  String gcs_output_vcf = gcs_output_dir + "/" + target_vcf
  String gcs_output_vcf_index = gcs_output_dir + "/" + target_vcf_index
  String gcs_output_sample_file = gcs_output_dir + "/" + target_sample_file

  command <<<

echo bcftools filter -R ~{include_region_bed} --threads ~{cpu} -o ~{target_vcf} ~{input_vcf}
bcftools filter -R ~{include_region_bed} --threads ~{cpu} -o ~{target_vcf} ~{input_vcf}

echo "build index"
bcftools index -t --threads ~{cpu} ~{target_vcf}

bcftools query -l ~{target_vcf} > ~{target_sample_file}

cat ~{target_sample_file} | wc -l > num_samples.txt

bcftools index -n ~{target_vcf} > num_variants.txt

echo gsutil ~{"-u " + project_id} -m cp ~{target_vcf} ~{gcs_output_vcf}
gsutil ~{"-u " + project_id} -m cp ~{target_vcf} ~{gcs_output_vcf}
rm -f ~{target_vcf}

echo gsutil ~{"-u " + project_id} -m cp ~{target_vcf_index} ~{gcs_output_vcf_index}
gsutil ~{"-u " + project_id} -m cp ~{target_vcf_index} ~{gcs_output_vcf_index}
rm -f ~{target_vcf_index}

echo gsutil ~{"-u " + project_id} -m cp ~{target_sample_file} ~{gcs_output_sample_file}
gsutil ~{"-u " + project_id} -m cp ~{target_sample_file} ~{gcs_output_sample_file}
rm -f ~{target_sample_file}

>>>

  runtime {
    docker: docker
    preemptible: preemptible
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    # The output has to be defined as String, otherwise the file would be delocalized and failed.
    String output_vcf = "~{gcs_output_vcf}"
    String output_vcf_index = "~{gcs_output_vcf_index}"
    String output_vcf_sample = "~{gcs_output_sample_file}"
    Int output_vcf_num_samples = read_int("num_samples.txt")
    Int output_vcf_num_variants = read_int("num_variants.txt")
  }
}
