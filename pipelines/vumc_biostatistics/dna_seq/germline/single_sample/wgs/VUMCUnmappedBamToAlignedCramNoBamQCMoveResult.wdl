version 1.0

import "../../../../../../tasks/broad/Utilities.wdl" as Utilities

workflow VUMCUnmappedBamToAlignedCramNoBamQCMoveResult {
  input {
    String genoset
    String GRID
    String target_bucket
    String? project_id

    Array[String] quality_yield_metrics

    String duplicate_metrics

    String? output_bqsr_reports

    String output_cram
    String output_cram_index
    String output_cram_md5
  }

  scatter(file in quality_yield_metrics){
    String moved_quality_yield_metrics_file = "~{target_bucket}/~{genoset}/~{GRID}/~{basename(file)}"
  }
  Array[String] moved_quality_yield_metrics = moved_quality_yield_metrics_file
  
  String moved_duplicate_metrics = "~{target_bucket}/~{genoset}/~{GRID}/~{basename(duplicate_metrics)}"

  String old_output_bqsr_reports = "~{output_bqsr_reports}"
  String moved_output_bqsr_reports = if old_output_bqsr_reports == "" then "" else "~{target_bucket}/~{genoset}/~{GRID}/~{basename(old_output_bqsr_reports)}"

  String moved_output_cram = "~{target_bucket}/~{genoset}/~{GRID}/~{basename(output_cram)}"
  String moved_output_cram_index = "~{target_bucket}/~{genoset}/~{GRID}/~{basename(output_cram_index)}"
  String moved_output_cram_md5 = "~{target_bucket}/~{genoset}/~{GRID}/~{basename(output_cram_md5)}"

  call MoveResult as mf {
    input:
      genoset = genoset,
      GRID = GRID,
      target_bucket = target_bucket,
      project_id = project_id,

      quality_yield_metrics = quality_yield_metrics,

      duplicate_metrics = duplicate_metrics,
      output_bqsr_reports = output_bqsr_reports,

      output_cram = output_cram,
      output_cram_index = output_cram_index,
      output_cram_md5 = output_cram_md5,

      moved_output_cram = moved_output_cram
  }

  output {
    Array[String] target_quality_yield_metrics = moved_quality_yield_metrics

    String target_duplicate_metrics = moved_duplicate_metrics

    String target_output_bqsr_reports = moved_output_bqsr_reports

    String target_output_cram = moved_output_cram
    String target_output_cram_index = moved_output_cram_index
    String target_output_cram_md5 = moved_output_cram_md5

    Int target_file_moved = mf.target_file_moved
  }
}

task MoveResult {
  input {
    String genoset
    String GRID
    String target_bucket
    String? project_id

    Array[String] quality_yield_metrics

    String duplicate_metrics
    String? output_bqsr_reports

    String output_cram
    String output_cram_index
    String output_cram_md5

    String moved_output_cram
  }
  

  command <<<
set +e

result=$(gsutil -q stat ~{output_cram} || echo 1)
if [[ $result != 1 ]]; then
  echo "Source cram file exists, moving to target bucket ..."

  set -e

  gsutil -m ~{"-u " + project_id} mv ~{sep=" " quality_yield_metrics} \
    ~{duplicate_metrics} \
    ~{output_bqsr_reports} \
    ~{output_cram} \
    ~{output_cram_index} \
    ~{output_cram_md5} \
    ~{target_bucket}/~{genoset}/~{GRID}/
else
  echo "Source cram file does not exist, checking target cram file ..."

  result=$(gsutil -q stat ~{moved_output_cram} || echo 1)
  if [[ $result != 1 ]]; then
    echo "Target cram file exists, return"
    exit 0
  else
    echo "Target cram file does not exist, error"
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
    Int target_file_moved = 1
  }
}
