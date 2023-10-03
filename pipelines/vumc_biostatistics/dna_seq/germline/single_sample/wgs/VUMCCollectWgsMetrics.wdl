version 1.0

import "../../../../../../tasks/broad/Qc.wdl" as QC
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"
import "./VUMCIlluminaCoverage.wdl" as VUMCIlluminaCoverage

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

  call CopyFileNoOverwrite as CopyFile {
    input:
      genoset = genoset,
      GRID = GRID,
      project_id = project_id,
      input_file = CollectWgsMetrics.metrics,
      target_bucket = target_bucket
  } 

  call VUMCIlluminaCoverage.GetIlluminaCoverage {
    input:
      wgs_metrics = CollectWgsMetrics.metrics
  }

  # Outputs that will be retained when execution is complete
  output {
    File wgs_metrics = CopyFile.output_file
    Float illumina_coverage = GetIlluminaCoverage.illumina_coverage
  }
  meta {
    allowNestedInputs: true
  }
}

task CopyFileNoOverwrite {
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

result=$(gsutil -q stat ~{new_file} || echo 1)
if [[ $result != 1 ]]; then
  echo "Target file exists, return"
  exit 0

  result=$(gsutil -q stat ~{input_file} || echo 1)
  if [[ $result != 1 ]]; then
    echo "Source file exists, copying to target bucket ..."

    set -e
      
    gsutil -m ~{"-u " + project_id} cp ~{input_file} \
      ~{target_bucket}/~{genoset}/~{GRID}/

  else
    echo "Both source file and target file does not exist, error ..."
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
