version 1.0

import "../genotype/Utils.wdl" as Utils

# For this workflow, we will expect the output is still plink bed format
workflow VUMCPlink2Filter {
  input {
    File source_bed
    File source_bim
    File source_fam

    String target_prefix

    String plink2_filter_option

    String? plink2_chr_option = "--chr-set 22 no-xy"

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

    String? project_id
    String? target_bucket
  }

  call Plink2Filter {
    input:
      source_bed = source_bed,
      source_bim = source_bim,
      source_fam = source_fam,

      plink2_filter_option = plink2_filter_option,
      plink2_chr_option = plink2_chr_option,

      target_prefix = target_prefix,

      docker = docker
  }

  if(defined(target_bucket)){
    call Utils.MoveOrCopyPlinkFile as CopyFile {
      input:
        source_bed = Plink2Filter.output_bed,
        source_bim = Plink2Filter.output_bim,
        source_fam = Plink2Filter.output_fam,
        is_move_file = false,
        project_id = project_id,
        target_bucket = select_first([target_bucket])
    }
  }

  output {
    File output_bed = select_first([CopyFile.output_bed, Plink2Filter.output_bed])
    File output_bim = select_first([CopyFile.output_bim, Plink2Filter.output_bim])
    File output_fam = select_first([CopyFile.output_fam, Plink2Filter.output_fam])
  }
}

task Plink2Filter {
  input {
      File source_bed
      File source_bim
      File source_fam

      String plink2_filter_option
      String? plink2_chr_option

      String target_prefix

      String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size(source_bed, "GB") * 2) + 2

  String new_bed = "${target_prefix}.bed"
  String new_bim = "${target_prefix}.bim"
  String new_fam = "${target_prefix}.fam"

  command <<<

plink2 \
  --bed ~{source_bed} \
  --bim ~{source_bim} \
  --fam ~{source_fam} \
  ~{plink2_chr_option} \
  ~{plink2_filter_option} \
  --make-bed \
  --out ~{target_prefix}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    File output_bed = new_bed
    File output_bim = new_bim
    File output_fam = new_fam
  }
}