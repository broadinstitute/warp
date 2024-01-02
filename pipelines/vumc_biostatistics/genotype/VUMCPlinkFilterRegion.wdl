version 1.0

import "./Utils.wdl" as Utils

workflow VUMCPlinkFilterRegion {
  input {
    File source_bed
    File source_bim
    File source_fam

    File region_bed
    String target_prefix
    String? target_suffix

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

    String? project_id
    String? target_bucket
  }

  String file_prefix=target_prefix + select_first([target_suffix, ""])

  call PlinkFilterRegion {
    input:
      source_bed = source_bed,
      source_bim = source_bim,
      source_fam = source_fam,
      region_bed = region_bed,
      target_prefix = file_prefix,
      docker = docker
  }

  if(defined(target_bucket)){
    call Utils.MoveOrCopyPlinkFile as CopyFile {
      input:
        source_bed = PlinkFilterRegion.output_bed,
        source_bim = PlinkFilterRegion.output_bim,
        source_fam = PlinkFilterRegion.output_fam,
        is_move_file = false,
        project_id = project_id,
        target_bucket = select_first([target_bucket])
    }
  }

  output {
    File output_bed = select_first([CopyFile.output_bed, PlinkFilterRegion.output_bed])
    File output_bim = select_first([CopyFile.output_bim, PlinkFilterRegion.output_bim])
    File output_fam = select_first([CopyFile.output_fam, PlinkFilterRegion.output_fam])
  }
}

task PlinkFilterRegion {
  input {
    File source_bed
    File source_bim
    File source_fam

    File region_bed
    
    String target_prefix
    
    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    Int memory_gb = 20
  }

  Int disk_size = ceil(size([source_bed, source_bim, source_fam], "GB")  * 2) + 20

  String new_bed = target_prefix + ".bed"
  String new_bim = target_prefix + ".bim"
  String new_fam = target_prefix + ".fam"

  command <<<

plink2 \
  --bed ~{source_bed} \
  --bim ~{source_bim} \
  --fam ~{source_fam} \
  --extract bed0 ~{region_bed} \
  --make-bed --out ~{target_prefix}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_bed = new_bed
    File output_bim = new_bim
    File output_fam = new_fam
  }
}
