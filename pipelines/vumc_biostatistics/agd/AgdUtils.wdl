version 1.0

task ReplaceICAIdWithGrid {
  input {
    File input_psam
    File id_map_file
    String output_psam
  }

  command <<<

python3 <<CODE

import io

with open("~{id_map_file}", "rt") as fin:
  id_map = {}
  for line in fin:
    parts = line.strip().split('\t')
    id_map[parts[0]] = parts[1]

with open("~{input_psam}", "rt") as fin:
  with open("~{output_psam}", "wt") as fout:
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
    File output_psam = "~{output_psam}"
  }
}

task CreateCohortPsam {
  input {
    File input_psam
    File grid_file
    String output_psam
  }

  command <<<

python3 <<CODE

import os

grids = set(line.strip() for line in open("~{grid_file}", "rt"))
with open("~{output_psam}", "wt") as fout:
  with open("~{input_psam}", "rt") as fin:
    for line in fin:
      if line.startswith("#"):
        fout.write(line)
        continue
      if line.split('\t')[1] in grids:
        fout.write(line)
CODE

echo "Number of samples to keep:"
wc -l "~{output_psam}"

>>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    File output_psam = "~{output_psam}"
  }
}
