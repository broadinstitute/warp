version 1.0

## Copyright Broad Institute/VUMC, 2018/2022
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in FASTQ format
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "../../../../../../tasks/broad/Qc.wdl" as QC
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow VUMCCollectWgsMetrics {

  String pipeline_version = "1.0.0.0"

  input {
    File input_cram
    File input_cram_index

    String? project_id
    String genoset
    String GRID
    String target_bucket

    DNASeqSingleSampleReferences references
    PapiSettings papi_settings

    File wgs_coverage_interval_list
  }

  # QC the sample WGS metrics (stringent thresholds)
  call QC.CollectWgsMetrics as CollectWgsMetrics {
    input:
      input_bam = input_cram,
      input_bam_index = input_cram_index,
      metrics_filename = GRID + ".wgs_metrics",
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  call MoveFile {
    input:
      genoset = genoset,
      GRID = GRID,
      project_id = project_id,
      input_file = CollectWgsMetrics.metrics,
      target_bucket = target_bucket
  } 

  call GetIlluminaCoverage {
    input:
      wgs_metrics = MoveFile.output_file
  }

  # Outputs that will be retained when execution is complete
  output {
    File wgs_metrics = MoveFile.output_file
    Float illumina_coverage = GetIlluminaCoverage.illumina_coverage
  }
  meta {
    allowNestedInputs: true
  }
}

task MoveFile {
  input {
    String genoset
    String GRID
    String? project_id

    String input_file

    String target_bucket
  }

  String new_file = "${target_bucket}/${genoset}/${GRID}/${basename(input_file)}"

  command <<<
set +e

result=$(gsutil -q stat ~{input_file} || echo 1)
if [[ $result != 1 ]]; then
  echo "Source file exists, moving to target bucket ..."

  set -e
    
  gsutil -m ~{"-u " + project_id} mv ~{input_file} \
    ~{target_bucket}/~{genoset}/~{GRID}/

else
  echo "Source file does not exist, checking target file ..."

  result=$(gsutil -q stat ~{new_file} || echo 1)
  if [[ $result != 1 ]]; then
    echo "Target file exists, return"
    exit 0
  else
    echo "Target file does not exist, error"
    exit 1
  fi
fi

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String output_file = "~{new_file}"
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
dat=fread(~{wgs_metrics})
tdat=t(dat)
coverage=tdat["MEAN_COVERAGE",1] * (1- tdat["PCT_EXC_DUPE",1]-tdat["PCT_EXC_OVERLAP",1])/(1-tdat["PCT_EXC_TOTAL",1])
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