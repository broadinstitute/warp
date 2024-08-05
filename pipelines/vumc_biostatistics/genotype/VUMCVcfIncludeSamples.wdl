version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils

workflow VUMCVcfIncludeSamples {
  input {
    File input_vcf
    File include_samples
    File? replace_samples
    String target_prefix
    String target_suffix = ".vcf.gz"
    String docker = "staphb/bcftools"

    String? project_id
    String? target_bucket
    String? genoset

    Int preemptible = 1
  }

  call BcftoolsIncludeSamples {
    input:
      input_vcf = input_vcf,
      include_samples = include_samples,
      replace_samples = replace_samples,
      target_prefix = target_prefix,
      target_suffix = target_suffix,
      docker = docker,
      preemptible = preemptible
  }

  if(defined(target_bucket)){
    String gcs_output_dir = sub("~{target_bucket}", "/+$", "")
    String target_gcp_folder = if(defined(genoset)) then "~{gcs_output_dir}/~{genoset}" else "~{gcs_output_dir}"

    call GcpUtils.MoveOrCopyThreeFiles {
      input:
        source_file1 = BcftoolsIncludeSamples.output_vcf,
        source_file2 = BcftoolsIncludeSamples.output_vcf_index,
        source_file3 = BcftoolsIncludeSamples.output_vcf_sample_file,
        is_move_file = true,
        project_id = project_id,
        target_gcp_folder = target_gcp_folder
    }
  }

  output {
    File output_vcf = select_first([MoveOrCopyThreeFiles.output_file1, BcftoolsIncludeSamples.output_vcf])
    File output_vcf_index = select_first([MoveOrCopyThreeFiles.output_file2, BcftoolsIncludeSamples.output_vcf_index])
    File output_vcf_sample_file = select_first([MoveOrCopyThreeFiles.output_file3, BcftoolsIncludeSamples.output_vcf_sample_file])
    Int num_samples = BcftoolsIncludeSamples.num_samples
    Int num_variants = BcftoolsIncludeSamples.num_variants
  }
}

task BcftoolsIncludeSamples {
  input {
    File input_vcf
    File include_samples
    File? replace_samples
    String target_prefix
    String target_suffix = ".vcf.gz"
    String docker = "staphb/bcftools"
    Float disk_factor = 3.0
    Int preemptible = 1
  }

  Int disk_size = ceil(size(input_vcf, "GB") * disk_factor) + 2
  String new_vcf = target_prefix + target_suffix
  Boolean has_replace_samples = "~{replace_samples}" != ""
  String output_sample_file = basename(input_vcf) + ".samples.txt"

  command <<<

bcftools head ~{input_vcf} > temp.vcf
zcat ~{input_vcf} | grep '^#CHROM' -m 1 -A 2 | grep -v '^#CHROM' >> temp.vcf

if [[ "~{has_replace_samples}" = "false" ]]; then
  echo "no replace_samples"
else
  #we need to get the header after replacing secondary grid with primary grid
  echo "has replace_samples, replace sample names in header"
  mv temp.vcf temp0.vcf
  bcftools reheader -s ~{replace_samples} temp0.vcf > temp.vcf
fi

echo "get all samples in original VCF"
bcftools head temp.vcf > header.txt
bcftools query -l temp.vcf > all.id.txt

echo "get included samples"
tr -d '\r' < ~{include_samples} > filter.id.txt
grep -Fxf all.id.txt filter.id.txt | sort | uniq > keep.id.txt

if [[ ! -s keep.id.txt ]]; then
  echo "ERROR: no samples to keep"
  exit 1
fi

zcat ~{input_vcf} | grep -v "^#" | cut -f 1-8 | head -n 10 > data.txt
if grep -Fq "Imputed" data.txt
then
  if grep -Fq "Imputed" header.txt
  then
    is_header_wrong=false
  else
    is_header_wrong=true
  fi
else
  is_header_wrong=false
fi

echo is_header_wrong = $is_header_wrong

if [ "$is_header_wrong" = "true" ];
then
  echo "header is wrong, need to add Imputed/Genotyped INFO in header"
  #there is error in merged imputation vcf files. The header has IMPUTED instead of Imputed.
  awk '/^#CHROM/ {printf("##INFO=<ID=Imputed,Number=0,Type=Flag,Description=\"Marker was imputed but NOT genotyped\">\n##INFO=<ID=Genotyped,Number=0,Type=Flag,Description=\"Marker was genotyped\">\n");} {print}' header.txt > new_header.txt
  bcftools reheader -h new_header.txt ~{input_vcf} | bcftools view -S keep.id.txt -o ~{new_vcf} -
else
  echo "header is correct"
  if [[ "~{has_replace_samples}" = "false" ]]; then
    echo "no replace_samples"
    echo bcftools view -S keep.id.txt -o ~{new_vcf} ~{input_vcf}
    bcftools view -S keep.id.txt -o ~{new_vcf} ~{input_vcf}
  else
    echo "has replace_samples"
    echo bcftools reheader -s ~{replace_samples} ~{input_vcf} | bcftools view -S keep.id.txt -o ~{new_vcf}
    bcftools reheader -s ~{replace_samples} ~{input_vcf} | bcftools view -S keep.id.txt -o ~{new_vcf}
  fi
fi

echo "build index"
bcftools index -t ~{new_vcf}

bcftools query -l ~{new_vcf} > ~{output_sample_file}
cat ~{output_sample_file} | wc -l > num_samples.txt

bcftools index -n ~{new_vcf} > num_variants.txt
>>>

  runtime {
    docker: docker
    preemptible: preemptible
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    File output_vcf = "~{new_vcf}"
    File output_vcf_index = "~{new_vcf}.tbi"
    File output_vcf_sample_file = "~{output_sample_file}"
    Int num_samples = read_int("num_samples.txt")
    Int num_variants = read_int("num_variants.txt")
  }
}
