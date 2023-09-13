version 1.0

workflow CombineGvsChromosome {

  String pipeline_version = "1.0.0"

  input {
    Array[File] intervals_files
    Array[String] vcf_files
  }

  call GetChromosomeVcfMap {
    input:
      intervals_files = intervals_files,
      vcf_files = vcf_files
  }

  scatter(chrom_vcf in GetChromosomeVcfMap.chrom_vcf_array) {
    call CombineGvsChromosome {
      input:
        chromosome = chrom_vcf.left,
        vcf_files = chrom_vcf.right
    }
  }

  output {
    Array[File] chromosome_vcfs = CombineGvsChromosome.output_vcf
  }
}

task GetChromosomeVcfMap {
  input {
    Array[File] intervals_files
    Array[String] vcf_files
    Int preemptible_tries=3
  }

  Int disk_size = ceil(size(intervals_files, "GB")) + 2

  command <<<

    python3 <<CODE

import pandas as pd
import os
import itertools
import json

from collections import OrderedDict

interval_files = "~{sep=',' intervals_files}".split(',')
vcf_files = "~{sep=',' vcf_files}".split(',')

df = [[pd.read_csv(f, sep='\t', header=None, comment='@').iloc[0][0], os.path.basename(f)] for f in interval_files]

key_func = lambda x: x[0]
chrom_map = OrderedDict()

for key, group in itertools.groupby(df, key_func):
  chrom_map[key] = [g[1].replace(".interval_list","") for g in group]

print(chrom_map.keys())

vcf_map = {os.path.basename(f):f for f in vcf_files}

result = [[chrom, [vcf_map[vcf_name] for vcf_name in chrom_map[chrom]]] for chrom in chrom_map]

with open('vcf_map.json', 'w') as fp:
  json.dump(result, fp)

CODE

exit 0

>>>

  runtime {
    docker: "pandas/pandas:mamba-minimal"
    preemptible: preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    Array[Pair[String, Array[String]]] chrom_vcf_array = read_json("vcf_map.json")
  }
}

task CombineGvsChromosome {
  input {
    String chromosome
    Array[File] vcf_files

    String docker = "staphb/bcftools"
    Int preemptible_tries=3
  }

  Int disk_size = ceil(size(vcf_files, "GB") * 2.5) + 2
  String chrom_vcf = chromosome + ".vcf.gz"
  String chrom_vcf_tbi = chromosome + ".vcf.gz.tbi"

  command <<<

bcftools concat -o ~{chrom_vcf} -O z ~{sep=' ' vcf_files}

tabix index ~{chrom_vcf}

exit 0

>>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    File output_vcf ="~{chrom_vcf}"
    File output_vcf_tbi = "~{chrom_vcf_tbi}"
  }
}