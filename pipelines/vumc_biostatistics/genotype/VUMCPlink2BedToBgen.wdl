version 1.0

import "./Utils.wdl" as Utils

workflow VUMCPlink2BedToBgen {
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

  call Plink2BedToBgen {
    input:
      source_bed = source_bed,
      source_bim = source_bim,
      source_fam = source_fam,
      target_prefix = target_prefix,

      docker = docker
  }

  if(defined(target_bucket)){
    call Utils.MoveOrCopyTwoFiles as CopyFile {
      input:
        source_file1 = Plink2BedToBgen.output_bgen,
        source_file2 = Plink2BedToBgen.output_sample,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_bucket])
    }
  }

  output {
    File output_bgen = select_first([CopyFile.output_file1, Plink2BedToBgen.output_bgen])
    File output_bgen_sample = select_first([CopyFile.output_file2, Plink2BedToBgen.output_sample])
  }
}

task Plink2BedToBgen {
  input {
    File source_bed
    File source_bim
    File source_fam
    
    String target_prefix
    
    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    Int memory_gb = 20
  }

  Int disk_size = ceil(size([source_bed, source_bim, source_fam], "GB")  * 2) + 20

  String new_bgen = target_prefix + ".bgen"
  String new_sample = target_prefix + ".sample"

  command <<<

## convert plink to pgen
plink2 \
  --bed ~{source_bed} \
  --bim ~{source_bim} \
  --fam ~{source_fam} \
  --export bgen-1.2 bits=8 \
  --out ~{target_prefix}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_bgen = new_bgen
    File output_sample = new_sample
  }
}
