version 1.0

workflow VUMCIlluminaCoverage {

  String pipeline_version = "1.0.0.0"

  input {
    File wgs_metrics
  }

  call GetIlluminaCoverage {
    input:
      wgs_metrics = wgs_metrics
  }

  output {
    Float illumina_coverage = GetIlluminaCoverage.illumina_coverage
  }
  meta {
    allowNestedInputs: true
  }
}

task GetIlluminaCoverage {
  input {
    File wgs_metrics
  }

  #https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-x-30x-coverage-technical-note-770-2014-042.pdf
  command <<<

R --vanilla <<RSCRIPT

library(data.table)
dat=fread("~{wgs_metrics}")
tdat=t(dat)
coverage=tdat["MEAN_COVERAGE",1] * (1 - tdat["PCT_EXC_DUPE",1] - tdat["PCT_EXC_OVERLAP",1]) / (1 - tdat["PCT_EXC_TOTAL",1])
writeLines(as.character(coverage), "coverage.txt")

RSCRIPT

>>>

  runtime {
    docker: "shengqh/cqs_scrnaseq:20230721"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    Float illumina_coverage = read_float("coverage.txt")
  }
}
