version 1.0

workflow VUMCVcfExtractSamples {
  input {
    File input_vcf
    File input_vcf_index
    File include_samples

    String target_prefix
    String target_suffix = ".vcf.gz"

    String? project_id
    String? target_gcp_folder
  }

  if(!defined(target_gcp_folder)){
    call BcftoolsExtractSamples {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        include_samples = include_samples,
        target_prefix = target_prefix,
        target_suffix = target_suffix
    }
  }

  if(defined(target_gcp_folder)){
    call BcftoolsExtractSamplesGcp {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        include_samples = include_samples,
        target_prefix = target_prefix,
        target_suffix = target_suffix,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_vcf = select_first([BcftoolsExtractSamples.output_vcf, BcftoolsExtractSamplesGcp.output_vcf])
    File output_vcf_index = select_first([BcftoolsExtractSamples.output_vcf_index, BcftoolsExtractSamplesGcp.output_vcf_index])
    File output_vcf_sample = select_first([BcftoolsExtractSamples.output_vcf_sample, BcftoolsExtractSamplesGcp.output_vcf_sample])
    Int output_vcf_num_samples = select_first([BcftoolsExtractSamples.output_vcf_num_samples, BcftoolsExtractSamplesGcp.output_vcf_num_samples])
    Int output_vcf_num_variants = select_first([BcftoolsExtractSamples.output_vcf_num_variants, BcftoolsExtractSamplesGcp.output_vcf_num_variants])
  }
}

task BcftoolsExtractSamples {
  input {
    File input_vcf
    File input_vcf_index
    File include_samples

    String target_prefix
    String target_suffix
    
    String docker = "shengqh/hail_gcp:20240213"
    Float disk_factor = 3.0
    Int preemptible = 1
  }

  Int disk_size = ceil(size(input_vcf, "GB") * disk_factor) + 2

  String target_vcf = target_prefix + target_suffix
  String target_vcf_index = target_vcf + ".tbi"
  String target_sample_file = target_vcf + ".samples.txt"

  command <<<

echo "get all samples in original VCF"
bcftools query -l ~{input_vcf} > all.id.txt

echo "get included samples"
tr -d '\r' < ~{include_samples} > filter.id.txt
grep -Fxf all.id.txt filter.id.txt | sort | uniq > keep.id.txt

if [[ ! -s keep.id.txt ]]; then
  echo "ERROR: no samples to keep"
  exit 1
fi

echo bcftools view -S keep.id.txt -o ~{target_vcf} ~{input_vcf}
bcftools view -S keep.id.txt -o ~{target_vcf} ~{input_vcf}

echo "build index"
bcftools index -t ~{target_vcf}

bcftools query -l ~{target_vcf} > ~{target_sample_file}

cat ~{target_sample_file} | wc -l > num_samples.txt

bcftools index -n ~{target_vcf} > num_variants.txt

>>>

  runtime {
    docker: docker
    preemptible: preemptible
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
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

task BcftoolsExtractSamplesGcp {
  input {
    File input_vcf
    File input_vcf_index
    File include_samples

    String target_prefix
    String target_suffix
    
    String? project_id
    String target_gcp_folder

    String docker = "shengqh/hail_gcp:20240213"
    Float disk_factor = 3.0
    Int preemptible = 1
  }

  Int disk_size = ceil(size(input_vcf, "GB") * disk_factor) + 2

  String target_vcf = target_prefix + target_suffix
  String target_vcf_index = target_vcf + ".tbi"
  String target_sample_file = target_vcf + ".samples.txt"

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")
  String gcs_output_vcf = gcs_output_dir + "/" + target_vcf
  String gcs_output_vcf_index = gcs_output_dir + "/" + target_vcf_index
  String gcs_output_sample_file = gcs_output_dir + "/" + target_sample_file

  command <<<

echo "get all samples in original VCF"
bcftools query -l ~{input_vcf} > all.id.txt

echo "get included samples"
tr -d '\r' < ~{include_samples} > filter.id.txt
grep -Fxf all.id.txt filter.id.txt | sort | uniq > keep.id.txt

if [[ ! -s keep.id.txt ]]; then
  echo "ERROR: no samples to keep"
  exit 1
fi

echo bcftools view -S keep.id.txt -o ~{target_vcf} ~{input_vcf}
bcftools view -S keep.id.txt -o ~{target_vcf} ~{input_vcf}

echo "build index"
bcftools index -t ~{target_vcf}

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
    memory: "2 GiB"
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
