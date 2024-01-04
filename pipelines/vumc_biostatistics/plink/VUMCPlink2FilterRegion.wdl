version 1.0

import "../genotype/Utils.wdl" as Utils

workflow VUMCPlink2FilterRegion {
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

  call Plink2FilterRegion {
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
        source_bed = Plink2FilterRegion.output_bed,
        source_bim = Plink2FilterRegion.output_bim,
        source_fam = Plink2FilterRegion.output_fam,
        is_move_file = false,
        project_id = project_id,
        target_bucket = select_first([target_bucket])
    }
  }

  File final_bed = select_first([CopyFile.output_bed, Plink2FilterRegion.output_bed])
  File final_bim = select_first([CopyFile.output_bim, Plink2FilterRegion.output_bim])
  File final_fam = select_first([CopyFile.output_fam, Plink2FilterRegion.output_fam])
  Float final_bed_size = size(final_bed)
  Float final_bim_size = size(final_bim)
  Float final_fam_size = size(final_fam)

  output {
    File output_bed = select_first([CopyFile.output_bed, Plink2FilterRegion.output_bed])
    File output_bim = select_first([CopyFile.output_bim, Plink2FilterRegion.output_bim])
    File output_fam = select_first([CopyFile.output_fam, Plink2FilterRegion.output_fam])
    Float output_bed_size = final_bed_size
    Float output_bim_size = final_bim_size
    Float output_fam_size = final_fam_size
  }
}

task Plink2FilterRegion {
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
