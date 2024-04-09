version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils

# Illumina VCF file index is csi format which might not be recognized by GATK tools
# We need to convert it to tbi format
workflow VUMCVcfIndex {
  input {
    File input_vcf

    String? project_id
    String? target_gcp_folder
  }

  call VcfIndex {
    input:
      input_vcf = input_vcf,
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyOneFile as CopyIndex {
      input:
        source_file = VcfIndex.output_vcf_index,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder]),
    }
  }

  output {
    File output_vcf_index = select_first([CopyIndex.output_file, VcfIndex.output_vcf_index])
  }
}

task VcfIndex {
  input {
    File input_vcf
    
    String docker = "shengqh/hail_gcp:20240213"

    Float disk_factor = 1.2
    Int preemptible = 1
    Int cpu = 1
  }

  Int disk_size = ceil(size(input_vcf, "GB") * disk_factor) + 10
  Int memory_gb = 2 * cpu

  String target_vcf_index = basename(input_vcf) + ".tbi"

  command <<<

echo "build index"
bcftools index -t --threads ~{cpu} ~{input_vcf}
cp ~{input_vcf}.tbi ~{target_vcf_index}

>>>

  runtime {
    docker: docker
    preemptible: preemptible
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    # The output has to be defined as File, otherwise the file would not be delocalized
    File output_vcf_index = "~{target_vcf_index}"
  }
}

