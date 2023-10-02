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

grep -A 1 GENOME_TERRITORY ~{wgs_metrics} > test.tsv

R --vanilla <<RSCRIPT

dat=read.table("test.tsv", header=T, sep="\t")
tdat=t(dat)

tn=as.numeric(tdat)
names(tn) = rownames(tdat)

coverage=tn["MEAN_COVERAGE"] * (1 - tn["PCT_EXC_DUPE"] - tn["PCT_EXC_OVERLAP"]) / (1 - tn["PCT_EXC_TOTAL"])
writeLines(as.character(coverage), "coverage.txt")

RSCRIPT

>>>

  runtime {
    docker: "r-base:4.3.1"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    Float illumina_coverage = read_float("coverage.txt")
  }
}
