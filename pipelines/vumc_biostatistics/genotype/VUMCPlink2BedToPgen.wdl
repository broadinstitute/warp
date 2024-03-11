version 1.0

import "./Utils.wdl" as Utils

workflow VUMCPlink2BedToPgen {
  input {
    File source_bed
    File source_bim
    File source_fam

    String target_prefix
    String? plink_option

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

    String? project_id
    String? target_bucket
  }

  call Plink2BedToPgen {
    input:
      source_bed = source_bed,
      source_bim = source_bim,
      source_fam = source_fam,
      target_prefix = target_prefix,

      docker = docker
  }

  if(defined(target_bucket)){
    call Utils.MoveOrCopyPlinkFile as CopyFile {
      input:
        source_bed = Plink2BedToPgen.output_pgen,
        source_bim = Plink2BedToPgen.output_pvar,
        source_fam = Plink2BedToPgen.output_psam,
        is_move_file = false,
        project_id = project_id,
        target_bucket = select_first([target_bucket])
    }
  }

  output {
    File output_pgen = select_first([CopyFile.output_bed, Plink2BedToPgen.output_pgen])
    File output_pvar = select_first([CopyFile.output_bim, Plink2BedToPgen.output_pvar])
    File output_psam = select_first([CopyFile.output_fam, Plink2BedToPgen.output_psam])
  }
}

task Plink2BedToPgen {
  input {
    File source_bed
    File source_bim
    File source_fam
    
    String target_prefix
    
    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    Int memory_gb = 20
  }

  Int disk_size = ceil(size([source_bed, source_bim, source_fam], "GB")  * 2) + 20

  String new_pgen = target_prefix + ".pgen"
  String new_pvar = target_prefix + ".pvar"
  String new_psam = target_prefix + ".psam"

  command <<<

## convert plink to pgen
plink2 \
  --bed ~{source_bed} \
  --bim ~{source_bim} \
  --fam ~{source_fam} \
  --make-pgen --out ~{target_prefix}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_pgen = new_pgen
    File output_pvar = new_pvar
    File output_psam = new_psam
  }
}
