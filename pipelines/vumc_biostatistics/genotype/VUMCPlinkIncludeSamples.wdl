version 1.0

import "./Utils.wdl" as Utils

workflow VUMCPlinkIncludeSamples {
  input {
    File source_bed
    File source_bim
    File source_fam

    File include_samples
    String target_prefix

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

    String? project_id
    String? target_bucket
  }

  call CreateIncludeFam {
    input:
      source_fam = source_fam,
      include_samples = include_samples,
  }

  call PlinkIncludeSamples {
    input:
      source_bed = source_bed,
      source_bim = source_bim,
      source_fam = source_fam,
      keep_id_fam = CreateIncludeFam.keep_id_fam,
      target_prefix = target_prefix,
      docker = docker
  }

  if(defined(target_bucket)){
    call Utils.MoveOrCopyPlinkFile {
      input:
        source_bed = PlinkIncludeSamples.output_bed,
        source_bim = PlinkIncludeSamples.output_bim,
        source_fam = PlinkIncludeSamples.output_fam,
        project_id = project_id,
        target_bucket = select_first([target_bucket]),
        action = "mv"
    }
  }

  output {
    File output_bed = select_first([MoveOrCopyPlinkFile.output_bed, PlinkIncludeSamples.output_bed])
    File output_bim = select_first([MoveOrCopyPlinkFile.output_bim, PlinkIncludeSamples.output_bim])
    File output_fam = select_first([MoveOrCopyPlinkFile.output_fam, PlinkIncludeSamples.output_fam])
  }
}

task CreateIncludeFam {
  input {
    File source_fam
    File include_samples
  }

  command <<<

python3 <<CODE

import os

grids = set(line.strip() for line in open("~{include_samples}", "rt"))
with open("keep.id.fam", "wt") as fout:
  with open("~{source_fam}", "rt") as fin:
    for line in fin:
      if line.split()[1] in grids:
        fout.write(line)
CODE

echo "Number of samples to keep:"
wc -l keep.id.fam

>>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    File keep_id_fam = "keep.id.fam"
  }
}

task PlinkIncludeSamples {
  input {
    File source_bed
    File source_bim
    File source_fam

    File keep_id_fam
    
    String target_prefix
    
    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([source_bed, source_bim, source_fam], "GB")  * 2) + 10

  String new_bed = target_prefix + ".bed"
  String new_bim = target_prefix + ".bim"
  String new_fam = target_prefix + ".fam"

  command <<<

plink2 \
  --bed ~{source_bed} \
  --bim ~{source_bim} \
  --fam ~{source_fam} \
  --keep ~{keep_id_fam} \
  --make-bed --out ~{target_prefix}

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
