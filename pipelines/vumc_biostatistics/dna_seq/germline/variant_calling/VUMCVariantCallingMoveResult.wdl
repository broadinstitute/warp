version 1.0

workflow VUMCVariantCallingMoveResult {
  input {
    String input_vcf_summary_metrics 
    String input_vcf_detail_metrics 
    String input_vcf
    String input_vcf_index 

    String? project_id
    String target_bucket
    String genoset
    String GRID
  }

  call CopyOrMoveResult as MoveResult {
    input: 
      input_vcf_summary_metrics = input_vcf_summary_metrics,
      input_vcf_detail_metrics = input_vcf_detail_metrics,
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      project_id = project_id,
      target_bucket = target_bucket,
      genoset = genoset,
      GRID = GRID,
      is_move_file = true
  }

  output {
    String vcf_summary_metrics = MoveResult.vcf_summary_metrics
    String vcf_detail_metrics = MoveResult.vcf_detail_metrics
    String output_vcf = MoveResult.output_vcf
    String output_vcf_index = MoveResult.output_vcf_index
  }
}

task CopyOrMoveResult {
  input {
    String input_vcf_summary_metrics 
    String input_vcf_detail_metrics 
    String input_vcf
    String input_vcf_index 

    Boolean is_move_file = true

    String? project_id
    String target_bucket
    String genoset
    String GRID
  }

  String action = if (is_move_file) then "mv" else "cp"

  String target_folder = "~{target_bucket}/~{genoset}/~{GRID}"

  String new_vcf_summary_metrics = "~{target_folder}/~{basename(input_vcf_summary_metrics)}"
  String new_vcf_detail_metrics = "~{target_folder}/~{basename(input_vcf_detail_metrics)}"
  String new_vcf = "~{target_folder}/~{basename(input_vcf)}"
  String new_vcf_index = "~{target_folder}/~{basename(input_vcf_index)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} \
  ~{input_vcf_summary_metrics} \
  ~{input_vcf_detail_metrics} \
  ~{input_vcf} \
  ~{input_vcf_index} \
  ~{target_folder}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String vcf_summary_metrics = new_vcf_summary_metrics
    String vcf_detail_metrics = new_vcf_detail_metrics
    String output_vcf = new_vcf
    String output_vcf_index = new_vcf_index
  }
}
