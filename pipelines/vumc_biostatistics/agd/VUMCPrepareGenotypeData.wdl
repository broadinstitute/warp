version 1.0

import "../genotype/Utils.wdl" as Utils

workflow VUMCPrepareGenotypeData {
  input {
    Array[File] source_pgen_files
    Array[File] source_psam_files
    Array[File] source_pvar_files

    Array[String] chromosomes

    String plink2_filter_option

    File grid_file
    String target_prefix

    File id_map_file

    String? project_id
    String? target_gcp_folder
  }

  scatter (idx in range(length(chromosomes))) {
    String chromosome = chromosomes[idx]
    File pgen_file = source_pgen_files[idx]
    File psam_file = source_psam_files[idx]
    File pvar_file = source_pvar_files[idx]
    String replaced_sample_name = "~{chromosome}.psam"

    call ReplaceIdSample {
      input:
        input_file = psam_file,
        id_map_file = id_map_file,
        output_filename = replaced_sample_name
    }

    String grid_sample_name = "~{chromosome}.grid.psam"
    call CreateCohortSample {
      input:
        input_file = ReplaceIdSample.output_file,
        grid_file = grid_file,
        output_filename = grid_sample_name
    }

    call PlinkExtractSamples {
      input:
        source_pgen = pgen_file,
        source_psam = ReplaceIdSample.output_file,
        source_pvar = pvar_file,
        chromosome = chromosome,
        plink2_filter_option = plink2_filter_option,
        extract_sample = CreateCohortSample.output_file
    }
  }

  call PlinkMergePgenFiles {
    input:
      pgen_files = PlinkExtractSamples.output_pgen_file,
      pvar_files = PlinkExtractSamples.output_pvar_file,
      psam_files = PlinkExtractSamples.output_psam_file,
      target_prefix = target_prefix
  }

  if(defined(target_gcp_folder)){
    call Utils.MoveOrCopyThreeFiles as CopyFile {
      input:
        source_file1 = PlinkMergePgenFiles.output_pgen_file,
        source_file2 = PlinkMergePgenFiles.output_pvar_file,
        source_file3 = PlinkMergePgenFiles.output_psam_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_pgen_file = select_first([CopyFile.output_file1, PlinkMergePgenFiles.output_pgen_file])
    File output_pvar_file = select_first([CopyFile.output_file2, PlinkMergePgenFiles.output_pvar_file])
    File output_psam_file = select_first([CopyFile.output_file3, PlinkMergePgenFiles.output_psam_file])
  }
}

task ReplaceIdSample {
  input {
    File input_file
    File id_map_file
    String output_filename
  }

  command <<<

python3 <<CODE

import gzip
import io

if "~{id_map_file}".endswith(".gz"):
  fin = gzip.open("~{id_map_file}", "rt")
else:
  fin = open("~{id_map_file}", "rt")

with fin:
  id_map = {}
  for line in fin:
    parts = line.strip().split('\t')
    id_map[parts[0]] = parts[1]

with open("~{input_file}", "rt") as fin:
  with open("~{output_filename}", "wt") as fout:
    for line in fin:
      parts = line.strip().split('\t')
      if parts[1] in id_map:
        grid = id_map[parts[1]]
        parts[1] = grid

      newline = '\t'.join(parts)
      fout.write(f"{newline}\n")

CODE

>>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    File output_file = "~{output_filename}"
  }
}

task CreateCohortSample {
  input {
    File input_file
    File grid_file
    String output_filename
  }

  command <<<

python3 <<CODE

import os

grids = set(line.strip() for line in open("~{grid_file}", "rt"))
with open("~{output_filename}", "wt") as fout:
  with open("~{input_file}", "rt") as fin:
    for line in fin:
      if line.startswith("#"):
        fout.write(line)
        continue
      if line.split('\t')[1] in grids:
        fout.write(line)
CODE

echo "Number of samples to keep:"
wc -l "~{output_filename}"

>>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    File output_file = "~{output_filename}"
  }
}

task PlinkExtractSamples {
  input {
    File source_pgen
    File source_psam
    File source_pvar
    File extract_sample
    String chromosome

    String plink2_filter_option

    Int memory_gb = 20

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([source_pgen, source_psam, source_pvar], "GB")  * 2) + 20

  String new_pgen = chromosome + ".pgen"
  String new_psam = chromosome + ".psam"
  String new_pvar = chromosome + ".pvar"

  command <<<

plink2 \
  --pgen ~{source_pgen} \
  --psam ~{source_psam} \
  --pvar ~{source_pvar} \
  ~{plink2_filter_option} \
  --keep ~{extract_sample} \
  --make-pgen \
  --out ~{chromosome}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_pgen_file = new_pgen
    File output_psam_file = new_psam
    File output_pvar_file = new_pvar
  }
}

task PlinkMergePgenFiles {
  input {
    Array[File] pgen_files
    Array[File] pvar_files
    Array[File] psam_files

    String target_prefix

    Int memory_gb = 20

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil((size(pgen_files, "GB") + size(pvar_files, "GB") + size(psam_files, "GB"))  * 3) + 20

  String new_pgen = target_prefix + ".pgen"
  String new_pvar = target_prefix + ".pvar"
  String new_psam = target_prefix + ".psam"

  String new_merged_pgen = target_prefix + "-merge.pgen"
  String new_merged_pvar = target_prefix + "-merge.pvar"
  String new_merged_psam = target_prefix + "-merge.psam"

  command <<<

cat ~{write_lines(pgen_files)} > pgen.list
cat ~{write_lines(pvar_files)} > pvar.list
cat ~{write_lines(psam_files)} > psam.list

paste pgen.list pvar.list psam.list > merge.list

plink2 --pmerge-list merge.list --make-pgen --out ~{target_prefix}

rm -f ~{new_pgen} ~{new_pvar} ~{new_psam}

mv ~{new_merged_pgen} ~{new_pgen}
mv ~{new_merged_pvar} ~{new_pvar}
mv ~{new_merged_psam} ~{new_psam}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_pgen_file = new_pgen
    File output_psam_file = new_psam
    File output_pvar_file = new_pvar
  }
}
