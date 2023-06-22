version 1.0

workflow VUMCMoveFastqToBamResult {
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

    String output_bam
    String output_bam_index
    String output_bam_md5

    String target_bucket
  }

  call MoveFastqToBamResult as mf {
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

      output_bam = output_bam,
      output_bam_index = output_bam_index,
      output_bam_md5 = output_bam_md5,

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

    String target_output_bam = mf.target_output_bam
    String target_output_bam_index = mf.target_output_bam_index
    String target_output_bam_md5 = mf.target_output_bam_md5

    Int target_bam_moved = mf.target_bam_moved
  }
}

task MoveFastqToBamResult {
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

    String output_bam
    String output_bam_index
    String output_bam_md5

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

  String new_output_bam = "${target_bucket}/${genoset}/${GRID}/${basename(output_bam)}"
  String new_output_bam_index = "${target_bucket}/${genoset}/${GRID}/${basename(output_bam_index)}"
  String new_output_bam_md5 = "${target_bucket}/${genoset}/${GRID}/${basename(output_bam_md5)}"

  command <<<

set -e

move_file(){
  SOURCE_FILE=$1
  TARGET_FILE=$2

  if [ $SOURCE_FILE != $TARGET_FILE ]; then
    gsutil -q stat $TARGET_FILE
    status=$?
    if [ $status -eq 0 ]; then
      echo "Target file exists, skipping move: $TARGET_FILE"
      return 0
    fi

    echo "Moving $SOURCE_FILE to $TARGET_FILE"
    gsutil mv $SOURCE_FILE $TARGET_FILE
    status=$?
    return $status    
  fi
}

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
move_file ~{output_bam} ~{new_output_bam}
move_file ~{output_bam_index} ~{new_output_bam_index}
move_file ~{output_bam_md5} ~{new_output_bam_md5}

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

    String target_output_bam = "~{new_output_bam}"
    String target_output_bam_index = "~{new_output_bam_index}"
    String target_output_bam_md5 = "~{new_output_bam_md5}"

    Int target_bam_moved = 1
  }
}
