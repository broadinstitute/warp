version 1.0

workflow TestOptionalValue {
  input {
    String duplicate_metrics
    String? output_bqsr_reports
  }

  call OptionalValue as mf {
    input:
      duplicate_metrics = duplicate_metrics,
      output_bqsr_reports = output_bqsr_reports,
  }

  output {
    String target_duplicate_metrics = mf.target_duplicate_metrics
    String? target_output_bqsr_reports = mf.target_output_bqsr_reports
  }
}

task OptionalValue {
  input {
    String duplicate_metrics
    String? output_bqsr_reports
  }

  String new_duplicate_metrics = "~{basename(duplicate_metrics)}"
  String bqsr = "~{output_bqsr_reports}"
  String new_output_bqsr_reports = if bqsr == "" then  "" else basename(bqsr)

  command <<<
move_file(){
  set +e

  echo "move_file $#"
  if [[ $# == 0 ]]; then
    echo "No arguments provided, skipping move."
    return 0
  fi
  
  SOURCE_FILE=$1
  TARGET_FILE=$2

  echo move $SOURCE_FILE $TARGET_FILE
  return 0
}
set -e

move_file ~{duplicate_metrics} ~{new_duplicate_metrics}

move_file ~{output_bqsr_reports} ~{new_output_bqsr_reports}

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_duplicate_metrics = "~{new_duplicate_metrics}"
    String? target_output_bqsr_reports = "~{new_output_bqsr_reports}"
  }
}
