version 1.0

import "./Utils.wdl" as Utils

workflow VUMCVcfExcludeSamples {
  input {
    File input_vcf
    File exclude_samples
    String target_prefix
    String target_suffix = ".vcf.gz"
    String docker = "staphb/bcftools"

    String? project_id
    String? target_bucket
    String? genoset
  }

  call BcftoolsExcludeSamples {
    input:
      input_vcf = input_vcf,
      exclude_samples = exclude_samples,
      target_prefix = target_prefix,
      target_suffix = target_suffix,
      docker = docker
  }

  if(defined(target_bucket)){
    call Utils.CopyVcfFile {
      input:
        input_vcf = BcftoolsExcludeSamples.output_vcf,
        input_vcf_index = BcftoolsExcludeSamples.output_vcf_index,
        project_id = project_id,
        target_bucket = select_first([target_bucket]),
        genoset = select_first([genoset]),
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
    File exclude_samples
    String target_prefix
    String target_suffix = ".vcf.gz"
    String docker = "staphb/bcftools"
  }

  Int disk_size = ceil(size(input_vcf, "GB") * 2) + 2
  String new_vcf = target_prefix + target_suffix

  command <<<

bcftools query -l ~{input_vcf} > all.id.txt

sort all.id.txt ~{exclude_samples} | uniq -d > remove.id.txt

bcftools head ~{input_vcf} > header.txt
zcat ~{input_vcf} | grep -v "^#" | cut -f 1-8 | head -n 1 > data.txt
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
  #there is error in merged imputation vcf files. The header has IMPUTED instead of Imputed.
  bcftools head ~{input_vcf} > header.txt
  awk '/^#CHROM/ {printf("##INFO=<ID=Imputed,Number=0,Type=Flag,Description=\"Marker was imputed but NOT genotyped\">\n##INFO=<ID=Genotyped,Number=0,Type=Flag,Description=\"Marker was genotyped\">\n");} {print}' header.txt > new_header.txt
  bcftools reheader -h new_header.txt ~{input_vcf} | bcftools view -S ^remove.id.txt -o ~{new_vcf} -
else
  bcftools view -S ^remove.id.txt -o ~{new_vcf} ~{input_vcf}
fi

bcftools index -t ~{new_vcf}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    File output_vcf = "~{new_vcf}"
    File output_vcf_index = "~{new_vcf}.tbi"
  }
}

