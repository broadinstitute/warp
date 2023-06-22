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
  if [ ~{unsorted_base_distribution_by_cycle_pdf} != ~{new_unsorted_base_distribution_by_cycle_pdf} ]; then
    gsutil mv ~{unsorted_base_distribution_by_cycle_pdf} ~{new_unsorted_base_distribution_by_cycle_pdf}
  fi
  
  if [ ~{unsorted_base_distribution_by_cycle_metrics} != ~{new_unsorted_base_distribution_by_cycle_metrics} ]; then
    gsutil mv ~{unsorted_base_distribution_by_cycle_metrics} ~{new_unsorted_base_distribution_by_cycle_metrics}
  fi
  
  if [ ~{unsorted_insert_size_histogram_pdf} != ~{new_unsorted_insert_size_histogram_pdf} ]; then
    gsutil mv ~{unsorted_insert_size_histogram_pdf} ~{new_unsorted_insert_size_histogram_pdf}
  fi

  if [ ~{unsorted_insert_size_metrics} != ~{new_unsorted_insert_size_metrics} ]; then
    gsutil mv ~{unsorted_insert_size_metrics} ~{new_unsorted_insert_size_metrics}
  fi

  if [ ~{unsorted_quality_by_cycle_pdf} != ~{new_unsorted_quality_by_cycle_pdf} ]; then
    gsutil mv ~{unsorted_quality_by_cycle_pdf} ~{new_unsorted_quality_by_cycle_pdf}
  fi

  if [ ~{unsorted_quality_by_cycle_metrics} != ~{new_unsorted_quality_by_cycle_metrics} ]; then
    gsutil mv ~{unsorted_quality_by_cycle_metrics} ~{new_unsorted_quality_by_cycle_metrics}
  fi

  if [ ~{unsorted_quality_distribution_pdf} != ~{new_unsorted_quality_distribution_pdf} ]; then
    gsutil mv ~{unsorted_quality_distribution_pdf} ~{new_unsorted_quality_distribution_pdf}
  fi

  if [ ~{unsorted_quality_distribution_metrics} != ~{new_unsorted_quality_distribution_metrics} ]; then
    gsutil mv ~{unsorted_quality_distribution_metrics} ~{new_unsorted_quality_distribution_metrics}
  fi

  if [ ~{cross_check_fingerprints_metrics} != ~{new_cross_check_fingerprints_metrics} ]; then
    gsutil mv ~{cross_check_fingerprints_metrics} ~{new_cross_check_fingerprints_metrics}
  fi

  if [ ~{selfSM} != ~{new_selfSM} ]; then
    gsutil mv ~{selfSM} ~{new_selfSM}
  fi

  if [ ~{duplicate_metrics} != ~{new_duplicate_metrics} ]; then
    gsutil mv ~{duplicate_metrics} ~{new_duplicate_metrics}
  fi
  
  if [ ~{output_bqsr_reports} != ~{new_output_bqsr_reports} ]; then
    gsutil mv ~{output_bqsr_reports} ~{new_output_bqsr_reports}
  fi

  if [ ~{output_cram} != ~{new_output_cram} ]; then
    gsutil mv ~{output_cram} ~{new_output_cram}
  fi

  if [ ~{output_cram_index} != ~{new_output_cram_index} ]; then
    gsutil mv ~{output_cram_index} ~{new_output_cram_index}
  fi

  if [ ~{output_cram_md5} != ~{new_output_cram_md5} ]; then
    gsutil mv ~{output_cram_md5} ~{new_output_cram_md5}
  fi
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
    String target_output_cram_md5 = "~{new_output_cram_index}"

    Int target_cram_moved = 1
  }
}
