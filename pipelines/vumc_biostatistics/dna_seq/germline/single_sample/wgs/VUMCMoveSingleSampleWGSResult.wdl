version 1.0

workflow VUMCMoveSingleSampleWGSResult {
  input {
    String genoset
    String GRID

    String read_group_alignment_summary_metrics
    String read_group_gc_bias_detail_metrics
    String read_group_gc_bias_pdf
    String read_group_gc_bias_summary_metrics

    String? cross_check_fingerprints_metrics

    String selfSM

    String agg_alignment_summary_metrics
    String agg_bait_bias_detail_metrics
    String agg_bait_bias_summary_metrics
    String agg_gc_bias_detail_metrics
    String agg_gc_bias_pdf
    String agg_gc_bias_summary_metrics
    String agg_insert_size_histogram_pdf
    String agg_insert_size_metrics
    String agg_pre_adapter_detail_metrics
    String agg_pre_adapter_summary_metrics
    String agg_quality_distribution_pdf
    String agg_quality_distribution_metrics
    String agg_error_summary_metrics

    String wgs_metrics
    String raw_wgs_metrics

    String duplicate_metrics
    String? output_bqsr_reports

    String output_cram
    String output_cram_index
    String output_cram_md5

    String validate_cram_file_report

    String output_vcf
    String output_vcf_index

    String gvcf_summary_metrics
    String gvcf_detail_metrics

    String target_bucket
  }

  call MoveSingleSampleWGSResult as mf {
    input:
      genoset = genoset,
      GRID = GRID,

      read_group_alignment_summary_metrics = read_group_alignment_summary_metrics,
      read_group_gc_bias_detail_metrics = read_group_gc_bias_detail_metrics,
      read_group_gc_bias_pdf = read_group_gc_bias_pdf,
      read_group_gc_bias_summary_metrics = read_group_gc_bias_summary_metrics,

      cross_check_fingerprints_metrics = cross_check_fingerprints_metrics,

      selfSM = selfSM,

      agg_alignment_summary_metrics = agg_alignment_summary_metrics,
      agg_bait_bias_detail_metrics = agg_bait_bias_detail_metrics,
      agg_bait_bias_summary_metrics = agg_bait_bias_summary_metrics,
      agg_gc_bias_detail_metrics = agg_gc_bias_detail_metrics,
      agg_gc_bias_pdf = agg_gc_bias_pdf,
      agg_gc_bias_summary_metrics = agg_gc_bias_summary_metrics,
      agg_insert_size_histogram_pdf = agg_insert_size_histogram_pdf,
      agg_insert_size_metrics = agg_insert_size_metrics,
      agg_pre_adapter_detail_metrics = agg_pre_adapter_detail_metrics,
      agg_pre_adapter_summary_metrics = agg_pre_adapter_summary_metrics,
      agg_quality_distribution_pdf = agg_quality_distribution_pdf,
      agg_quality_distribution_metrics = agg_quality_distribution_metrics,
      agg_error_summary_metrics = agg_error_summary_metrics,

      wgs_metrics = wgs_metrics,
      raw_wgs_metrics = raw_wgs_metrics,

      duplicate_metrics = duplicate_metrics,
      output_bqsr_reports = output_bqsr_reports,

      output_cram = output_cram,
      output_cram_index = output_cram_index,
      output_cram_md5 = output_cram_md5,

      validate_cram_file_report = validate_cram_file_report,

      output_vcf = output_vcf,
      output_vcf_index = output_vcf_index,

      gvcf_summary_metrics = gvcf_summary_metrics,
      gvcf_detail_metrics = gvcf_detail_metrics,

      target_bucket = target_bucket
  }

  output {
    String target_read_group_alignment_summary_metrics = mf.target_read_group_alignment_summary_metrics
    String target_read_group_gc_bias_detail_metrics = mf.target_read_group_gc_bias_detail_metrics
    String target_read_group_gc_bias_pdf = mf.target_read_group_gc_bias_pdf
    String target_read_group_gc_bias_summary_metrics = mf.target_read_group_gc_bias_summary_metrics

    String? target_cross_check_fingerprints_metrics = mf.target_cross_check_fingerprints_metrics

    String target_selfSM = mf.target_selfSM

    String target_agg_alignment_summary_metrics = mf.target_agg_alignment_summary_metrics
    String target_agg_bait_bias_detail_metrics = mf.target_agg_bait_bias_detail_metrics
    String target_agg_bait_bias_summary_metrics = mf.target_agg_bait_bias_summary_metrics
    String target_agg_gc_bias_detail_metrics = mf.target_agg_gc_bias_detail_metrics
    String target_agg_gc_bias_pdf = mf.target_agg_gc_bias_pdf
    String target_agg_gc_bias_summary_metrics = mf.target_agg_gc_bias_summary_metrics
    String target_agg_insert_size_histogram_pdf = mf.target_agg_insert_size_histogram_pdf
    String target_agg_insert_size_metrics = mf.target_agg_insert_size_metrics
    String target_agg_pre_adapter_detail_metrics = mf.target_agg_pre_adapter_detail_metrics
    String target_agg_pre_adapter_summary_metrics = mf.target_agg_pre_adapter_summary_metrics
    String target_agg_quality_distribution_pdf = mf.target_agg_quality_distribution_pdf
    String target_agg_quality_distribution_metrics = mf.target_agg_quality_distribution_metrics
    String target_agg_error_summary_metrics = mf.target_agg_error_summary_metrics

    String target_wgs_metrics = mf.target_wgs_metrics
    String target_raw_wgs_metrics = mf.target_raw_wgs_metrics

    String target_duplicate_metrics = mf.target_duplicate_metrics
    String? target_output_bqsr_reports = mf.target_output_bqsr_reports

    String target_output_cram = mf.target_output_cram
    String target_output_cram_index = mf.target_output_cram_index
    String target_output_cram_md5 = mf.target_output_cram_md5

    String target_validate_cram_file_report = mf.target_validate_cram_file_report

    String target_output_vcf = mf.target_output_vcf
    String target_output_vcf_index = mf.target_output_vcf_index

    String target_gvcf_summary_metrics = mf.target_gvcf_summary_metrics
    String target_gvcf_detail_metrics = mf.target_gvcf_detail_metrics

    Int target_file_moved = mf.target_file_moved
  }
}

task MoveSingleSampleWGSResult {
  input {
    String genoset
    String GRID

    String read_group_alignment_summary_metrics
    String read_group_gc_bias_detail_metrics
    String read_group_gc_bias_pdf
    String read_group_gc_bias_summary_metrics

    String? cross_check_fingerprints_metrics

    String selfSM

    String agg_alignment_summary_metrics
    String agg_bait_bias_detail_metrics
    String agg_bait_bias_summary_metrics
    String agg_gc_bias_detail_metrics
    String agg_gc_bias_pdf
    String agg_gc_bias_summary_metrics
    String agg_insert_size_histogram_pdf
    String agg_insert_size_metrics
    String agg_pre_adapter_detail_metrics
    String agg_pre_adapter_summary_metrics
    String agg_quality_distribution_pdf
    String agg_quality_distribution_metrics
    String agg_error_summary_metrics

    String wgs_metrics
    String raw_wgs_metrics

    String duplicate_metrics
    String? output_bqsr_reports

    String output_cram
    String output_cram_index
    String output_cram_md5

    String validate_cram_file_report

    String output_vcf
    String output_vcf_index

    String gvcf_summary_metrics
    String gvcf_detail_metrics

    String target_bucket
  }

  String new_read_group_alignment_summary_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(read_group_alignment_summary_metrics)}"
  String new_read_group_gc_bias_detail_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(read_group_gc_bias_detail_metrics)}"
  String new_read_group_gc_bias_pdf = "${target_bucket}/${genoset}/${GRID}/${basename(read_group_gc_bias_pdf)}"
  String new_read_group_gc_bias_summary_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(read_group_gc_bias_summary_metrics)}"

  String? new_cross_check_fingerprints_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(cross_check_fingerprints_metrics)}"

  String? move_cross_check_fingerprints_metrics = "move_file ${cross_check_fingerprints_metrics} ${new_cross_check_fingerprints_metrics}"

  String new_selfSM = "${target_bucket}/${genoset}/${GRID}/${basename(selfSM)}"

  String new_agg_alignment_summary_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_alignment_summary_metrics)}"
  String new_agg_bait_bias_detail_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_bait_bias_detail_metrics)}"
  String new_agg_bait_bias_summary_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_bait_bias_summary_metrics)}"
  String new_agg_gc_bias_detail_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_gc_bias_detail_metrics)}"
  String new_agg_gc_bias_pdf = "${target_bucket}/${genoset}/${GRID}/${basename(agg_gc_bias_pdf)}"
  String new_agg_gc_bias_summary_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_gc_bias_summary_metrics)}"
  String new_agg_insert_size_histogram_pdf = "${target_bucket}/${genoset}/${GRID}/${basename(agg_insert_size_histogram_pdf)}"
  String new_agg_insert_size_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_insert_size_metrics)}"
  String new_agg_pre_adapter_detail_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_pre_adapter_detail_metrics)}"
  String new_agg_pre_adapter_summary_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_pre_adapter_summary_metrics)}"
  String new_agg_quality_distribution_pdf = "${target_bucket}/${genoset}/${GRID}/${basename(agg_quality_distribution_pdf)}"
  String new_agg_quality_distribution_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_quality_distribution_metrics)}"
  String new_agg_error_summary_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(agg_error_summary_metrics)}"

  String new_wgs_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(wgs_metrics)}"
  String new_raw_wgs_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(raw_wgs_metrics)}"

  String new_duplicate_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(duplicate_metrics)}"
  String? new_output_bqsr_reports = "${target_bucket}/${genoset}/${GRID}/${basename(output_bqsr_reports)}"

  String? move_output_bqsr_reports = "move_file ${output_bqsr_reports} ${new_output_bqsr_reports}"

  String new_output_cram = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram)}"
  String new_output_cram_index = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram_index)}"
  String new_output_cram_md5 = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram_md5)}"

  String new_validate_cram_file_report = "${target_bucket}/${genoset}/${GRID}/${basename(validate_cram_file_report)}"

  String new_output_vcf = "${target_bucket}/${genoset}/${GRID}/${basename(output_vcf)}"
  String new_output_vcf_index = "${target_bucket}/${genoset}/${GRID}/${basename(output_vcf_index)}"
  
  String new_gvcf_summary_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(gvcf_summary_metrics)}"
  String new_gvcf_detail_metrics = "${target_bucket}/${genoset}/${GRID}/${basename(gvcf_detail_metrics)}"

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

move_file ~{read_group_alignment_summary_metrics} ~{new_read_group_alignment_summary_metrics}
move_file ~{read_group_gc_bias_detail_metrics} ~{new_read_group_gc_bias_detail_metrics}
move_file ~{read_group_gc_bias_pdf} ~{new_read_group_gc_bias_pdf}
move_file ~{read_group_gc_bias_summary_metrics} ~{new_read_group_gc_bias_summary_metrics}

~{move_cross_check_fingerprints_metrics}

move_file ~{selfSM} ~{new_selfSM}

move_file ~{agg_alignment_summary_metrics} ~{new_agg_alignment_summary_metrics}
move_file ~{agg_bait_bias_detail_metrics} ~{new_agg_bait_bias_detail_metrics}
move_file ~{agg_bait_bias_summary_metrics} ~{new_agg_bait_bias_summary_metrics}
move_file ~{agg_gc_bias_detail_metrics} ~{new_agg_gc_bias_detail_metrics}
move_file ~{agg_gc_bias_pdf} ~{new_agg_gc_bias_pdf}
move_file ~{agg_gc_bias_summary_metrics} ~{new_agg_gc_bias_summary_metrics}
move_file ~{agg_insert_size_histogram_pdf} ~{new_agg_insert_size_histogram_pdf}
move_file ~{agg_insert_size_metrics} ~{new_agg_insert_size_metrics}
move_file ~{agg_pre_adapter_detail_metrics} ~{new_agg_pre_adapter_detail_metrics}
move_file ~{agg_pre_adapter_summary_metrics} ~{new_agg_pre_adapter_summary_metrics}
move_file ~{agg_quality_distribution_pdf} ~{new_agg_quality_distribution_pdf}
move_file ~{agg_quality_distribution_metrics} ~{new_agg_quality_distribution_metrics}
move_file ~{agg_error_summary_metrics} ~{new_agg_error_summary_metrics}

move_file ~{wgs_metrics} ~{new_wgs_metrics}
move_file ~{raw_wgs_metrics} ~{new_raw_wgs_metrics}

move_file ~{duplicate_metrics} ~{new_duplicate_metrics}

~{move_output_bqsr_reports}

move_file ~{output_cram} ~{new_output_cram}
move_file ~{output_cram_index} ~{new_output_cram_index}
move_file ~{output_cram_md5} ~{new_output_cram_md5}

move_file ~{validate_cram_file_report} ~{new_validate_cram_file_report}

move_file ~{output_vcf} ~{new_output_vcf}
move_file ~{output_vcf_index} ~{new_output_vcf_index}

move_file ~{gvcf_summary_metrics} ~{new_gvcf_summary_metrics}
move_file ~{gvcf_detail_metrics} ~{new_gvcf_detail_metrics}

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_read_group_alignment_summary_metrics = "~{new_read_group_alignment_summary_metrics}"
    String target_read_group_gc_bias_detail_metrics = "~{new_read_group_gc_bias_detail_metrics}"
    String target_read_group_gc_bias_pdf = "~{new_read_group_gc_bias_pdf}"
    String target_read_group_gc_bias_summary_metrics = "~{new_read_group_gc_bias_summary_metrics}"

    String? target_cross_check_fingerprints_metrics = "~{new_cross_check_fingerprints_metrics}"

    String target_selfSM = "~{new_selfSM}"

    String target_agg_alignment_summary_metrics = "~{new_agg_alignment_summary_metrics}"
    String target_agg_bait_bias_detail_metrics = "~{new_agg_bait_bias_detail_metrics}"
    String target_agg_bait_bias_summary_metrics = "~{new_agg_bait_bias_summary_metrics}"
    String target_agg_gc_bias_detail_metrics = "~{new_agg_gc_bias_detail_metrics}"
    String target_agg_gc_bias_pdf = "~{new_agg_gc_bias_pdf}"
    String target_agg_gc_bias_summary_metrics = "~{new_agg_gc_bias_summary_metrics}"
    String target_agg_insert_size_histogram_pdf = "~{new_agg_insert_size_histogram_pdf}"
    String target_agg_insert_size_metrics = "~{new_agg_insert_size_metrics}"
    String target_agg_pre_adapter_detail_metrics = "~{new_agg_pre_adapter_detail_metrics}"
    String target_agg_pre_adapter_summary_metrics = "~{new_agg_pre_adapter_summary_metrics}"
    String target_agg_quality_distribution_pdf = "~{new_agg_quality_distribution_pdf}"
    String target_agg_quality_distribution_metrics = "~{new_agg_quality_distribution_metrics}"
    String target_agg_error_summary_metrics = "~{new_agg_error_summary_metrics}"

    String target_wgs_metrics = "~{new_wgs_metrics}"
    String target_raw_wgs_metrics = "~{new_raw_wgs_metrics}"

    String target_duplicate_metrics = "~{new_duplicate_metrics}"
    String? target_output_bqsr_reports = "~{new_output_bqsr_reports}"

    String target_output_cram = "~{new_output_cram}"
    String target_output_cram_index = "~{new_output_cram_index}"
    String target_output_cram_md5 = "~{new_output_cram_md5}"

    String target_validate_cram_file_report = "~{new_validate_cram_file_report}"

    String target_output_vcf = "~{new_output_vcf}"
    String target_output_vcf_index = "~{new_output_vcf_index}"

    String target_gvcf_summary_metrics = "~{new_gvcf_summary_metrics}"
    String target_gvcf_detail_metrics = "~{new_gvcf_detail_metrics}"

    Int target_file_moved = 1
  }
}
