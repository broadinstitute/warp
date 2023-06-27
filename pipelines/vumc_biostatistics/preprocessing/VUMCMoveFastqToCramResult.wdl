version 1.0

workflow VUMCMoveFastqToCramResult {
  input {
    String genoset
    String GRID

    String unsorted_base_distribution_by_cycle_pdf
    String unsorted_base_distribution_by_cycle_metrics
    String unsorted_insert_size_histogram_pdf
    String unsorted_insert_size_metrics
    String unsorted_quality_by_cycle_pdf
    String unsorted_quality_by_cycle_metrics
    String unsorted_quality_distribution_pdf
    String unsorted_quality_distribution_metrics

    String cross_check_fingerprints_metrics

    String selfSM

    String duplicate_metrics
    String output_bqsr_reports

    String output_cram
    String output_cram_index
    String output_cram_md5

    String target_bucket
  }

  call MoveFastqToCramResult as mf {
    input:
      genoset = genoset,
      GRID = GRID,

      unsorted_base_distribution_by_cycle_pdf = unsorted_base_distribution_by_cycle_pdf, 
      unsorted_base_distribution_by_cycle_metrics = unsorted_base_distribution_by_cycle_metrics,
      unsorted_insert_size_histogram_pdf = unsorted_insert_size_histogram_pdf,
      unsorted_insert_size_metrics = unsorted_insert_size_metrics,
      unsorted_quality_by_cycle_pdf = unsorted_quality_by_cycle_pdf,
      unsorted_quality_by_cycle_metrics = unsorted_quality_by_cycle_metrics,
      unsorted_quality_distribution_pdf = unsorted_quality_distribution_pdf,
      unsorted_quality_distribution_metrics = unsorted_quality_distribution_metrics,

      cross_check_fingerprints_metrics = cross_check_fingerprints_metrics,

      selfSM = selfSM,

      duplicate_metrics = duplicate_metrics,
      output_bqsr_reports = output_bqsr_reports,

      output_cram = output_cram,
      output_cram_index = output_cram_index,
      output_cram_md5 = output_cram_md5,

      target_bucket = target_bucket
  }

  output {
    String target_unsorted_base_distribution_by_cycle_pdf = mf.target_unsorted_base_distribution_by_cycle_pdf
    String target_unsorted_base_distribution_by_cycle_metrics = mf.target_unsorted_base_distribution_by_cycle_metrics
    String target_unsorted_insert_size_histogram_pdf = mf.target_unsorted_insert_size_histogram_pdf
    String target_unsorted_insert_size_metrics = mf.target_unsorted_insert_size_metrics
    String target_unsorted_quality_by_cycle_pdf = mf.target_unsorted_quality_by_cycle_pdf
    String target_unsorted_quality_by_cycle_metrics = mf.target_unsorted_quality_by_cycle_metrics
    String target_unsorted_quality_distribution_pdf = mf.target_unsorted_quality_distribution_pdf
    String target_unsorted_quality_distribution_metrics = mf.target_unsorted_quality_distribution_metrics

    String target_cross_check_fingerprints_metrics = mf.target_cross_check_fingerprints_metrics

    String target_selfSM = mf.target_selfSM

    String target_duplicate_metrics = mf.target_duplicate_metrics
    String target_output_bqsr_reports = mf.target_output_bqsr_reports

    String target_output_cram = mf.target_output_cram
    String target_output_cram_index = mf.target_output_cram_index
    String target_output_cram_md5 = mf.target_output_cram_md5

    Int target_cram_moved = mf.target_cram_moved
  }
}

task MoveFastqToCramResult {
  input {
    String genoset
    String GRID

    String unsorted_base_distribution_by_cycle_pdf
    String unsorted_base_distribution_by_cycle_metrics
    String unsorted_insert_size_histogram_pdf
    String unsorted_insert_size_metrics
    String unsorted_quality_by_cycle_pdf
    String unsorted_quality_by_cycle_metrics
    String unsorted_quality_distribution_pdf
    String unsorted_quality_distribution_metrics

    String cross_check_fingerprints_metrics

    String selfSM

    String duplicate_metrics
    String output_bqsr_reports

    String output_cram
    String output_cram_index
    String output_cram_md5

    String target_bucket
  }

  String new_unsorted_base_distribution_by_cycle_pdf = "${target_bucket}/${genoset}/${GRID}/${basename(unsorted_base_distribution_by_cycle_pdf)}"
  String new_unsorted_base_distribution_by_cycle_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(unsorted_base_distribution_by_cycle_metrics)}"
  String new_unsorted_insert_size_histogram_pdf = "${target_bucket}/${genoset}/${GRID}/${basename(unsorted_insert_size_histogram_pdf)}"
  String new_unsorted_insert_size_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(unsorted_insert_size_metrics)}"
  String new_unsorted_quality_by_cycle_pdf = "${target_bucket}/${genoset}/${GRID}/${basename(unsorted_quality_by_cycle_pdf)}"
  String new_unsorted_quality_by_cycle_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(unsorted_quality_by_cycle_metrics)}"
  String new_unsorted_quality_distribution_pdf = "${target_bucket}/${genoset}/${GRID}/${basename(unsorted_quality_distribution_pdf)}"
  String new_unsorted_quality_distribution_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(unsorted_quality_distribution_metrics)}"

  String new_cross_check_fingerprints_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(cross_check_fingerprints_metrics)}"

  String new_selfSM = "${target_bucket}/${genoset}/${GRID}/${basename(selfSM)}"

  String new_duplicate_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(duplicate_metrics)}"
  String new_output_bqsr_reports = "${target_bucket}/${genoset}/${GRID}/${basename(output_bqsr_reports)}"

  String new_output_cram = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram)}"
  String new_output_cram_index = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram_index)}"
  String new_output_cram_md5 = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram_md5)}"

  command <<<
  

move_file(){
  set +e

  SOURCE_FILE=$1
  TARGET_FILE=$2

  echo "Moving $SOURCE_FILE to $TARGET_FILE"

  if [[ $SOURCE_FILE == $TARGET_FILE ]]; then
    echo "Target file equals to source file, skipping move: $TARGET_FILE"
    return 0
  fi

  echo "Checking if target file exists: $TARGET_FILE"

  gsutil -q stat $TARGET_FILE
  status=$?
  if [[ $status -eq 0 ]]; then
    echo "Target file exists, skipping move: $TARGET_FILE"
    return 0
  fi

  echo gsutil mv $SOURCE_FILE $TARGET_FILE

  gsutil mv $SOURCE_FILE $TARGET_FILE
  status=$?
  if [[ $status -eq 0 ]]; then
    echo "Moving succeed."
  else
    echo "Moving failed with status: $status"
  fi

  set -e
  return $status
}

set -e

move_file ~{unsorted_base_distribution_by_cycle_pdf} ~{new_unsorted_base_distribution_by_cycle_pdf} 
move_file ~{unsorted_base_distribution_by_cycle_metrics} ~{new_unsorted_base_distribution_by_cycle_metrics} 
move_file ~{unsorted_insert_size_histogram_pdf} ~{new_unsorted_insert_size_histogram_pdf}
move_file ~{unsorted_insert_size_metrics} ~{new_unsorted_insert_size_metrics}
move_file ~{unsorted_quality_by_cycle_pdf} ~{new_unsorted_quality_by_cycle_pdf}
move_file ~{unsorted_quality_by_cycle_metrics} ~{new_unsorted_quality_by_cycle_metrics}
move_file ~{unsorted_quality_distribution_pdf} ~{new_unsorted_quality_distribution_pdf}
move_file ~{unsorted_quality_distribution_metrics} ~{new_unsorted_quality_distribution_metrics}
move_file ~{cross_check_fingerprints_metrics} ~{new_cross_check_fingerprints_metrics}
move_file ~{selfSM} ~{new_selfSM}
move_file ~{duplicate_metrics} ~{new_duplicate_metrics}
move_file ~{output_bqsr_reports} ~{new_output_bqsr_reports}
move_file ~{output_cram} ~{new_output_cram}
move_file ~{output_cram_index} ~{new_output_cram_index}
move_file ~{output_cram_md5} ~{new_output_cram_md5}

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_unsorted_base_distribution_by_cycle_pdf = "~{new_unsorted_base_distribution_by_cycle_pdf}"
    String target_unsorted_base_distribution_by_cycle_metrics = "~{new_unsorted_base_distribution_by_cycle_metrics}"
    String target_unsorted_insert_size_histogram_pdf = "~{new_unsorted_insert_size_histogram_pdf}"
    String target_unsorted_insert_size_metrics = "~{new_unsorted_insert_size_metrics}"
    String target_unsorted_quality_by_cycle_pdf = "~{new_unsorted_quality_by_cycle_pdf}"
    String target_unsorted_quality_by_cycle_metrics = "~{new_unsorted_quality_by_cycle_metrics}"
    String target_unsorted_quality_distribution_pdf = "~{new_unsorted_quality_distribution_pdf}"
    String target_unsorted_quality_distribution_metrics = "~{new_unsorted_quality_distribution_metrics}"

    String target_cross_check_fingerprints_metrics = "~{new_cross_check_fingerprints_metrics}"

    String target_selfSM = "~{new_selfSM}"

    String target_duplicate_metrics = "~{new_duplicate_metrics}"
    String target_output_bqsr_reports = "~{new_output_bqsr_reports}"

    String target_output_cram = "~{new_output_cram}"
    String target_output_cram_index = "~{new_output_cram_index}"
    String target_output_cram_md5 = "~{new_output_cram_md5}"

    Int target_cram_moved = 1
  }
}
