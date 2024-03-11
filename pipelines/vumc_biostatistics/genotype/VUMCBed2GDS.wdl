version 1.0

import "./Utils.wdl" as Utils

workflow VUMCBed2GDS {
  input {
    File source_bed
    File source_bim
    File source_fam

    String reference_genome = "GRCh38"

    String target_prefix

    String docker = "shengqh/gds:20240205"

    Int memory_gb = 20

    String? project_id
    String? target_gcp_folder
  }

  call Bed2GDS {
    input:
      source_bed = source_bed,
      source_bim = source_bim,
      source_fam = source_fam,
      target_prefix = target_prefix,
      docker = docker,
      memory_gb = memory_gb
  }

  if(defined(target_gcp_folder)){
    call Utils.MoveOrCopyOneFile as CopyFile {
      input:
        source_file = Bed2GDS.gds_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File gds_file = select_first([CopyFile.output_file, Bed2GDS.gds_file])
  }
}

task Bed2GDS {
  input {
    File source_bed
    File source_bim
    File source_fam
    String target_prefix
    String docker
    Int memory_gb
  }

  Int disk_size = ceil(size([source_bed, source_bim, source_fam], "GB")  * 1.5) + 20
  Int total_memory_gb = memory_gb + 2
  String target_file = "~{target_prefix}.gds"

  command <<<

cat <<EOT >> convertBEDtoGDS.r
#############################################################################
#Title: convertBEDtoGDS
#Function:
# * Build the GDS file from plink binary files
#Author: Jiang, Lan/Sheng, Quanhu (modified from convertVCFtoGDS.r
#Time:     
#############################################################################
library(gdsfmt)
library(SeqArray)

seqBED2GDS("~{source_bed}", "~{source_fam}", "~{source_bim}", "~{target_file}", compress.geno="LZMA_RA",compress.annotation="LZMA_RA")
print("GDS build")

genofile<-seqOpen("~{target_file}")
print("GDS details")

###Closing Up###
genofile
seqClose(genofile)

EOT

R --vanilla -f convertBEDtoGDS.r

>>>

  runtime {
    docker: "~{docker}"
    preemptible: 1
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{total_memory_gb} GiB"
  }
  output {
    File gds_file = "~{target_file}"
  }
}
