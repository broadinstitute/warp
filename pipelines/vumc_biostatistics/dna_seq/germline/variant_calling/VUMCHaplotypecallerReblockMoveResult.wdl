version 1.0

workflow VUMCHaplotypecallerReblockMoveResult {
  input {
    String genoset
    String GRID
    String? project_id

    String input_vcf
    String input_vcf_index

    String target_bucket
  }

  call MoveVcf as mh {
    input:
      genoset = genoset,
      GRID = GRID,
      project_id = project_id,

      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,

      target_bucket = target_bucket
  }

  output {
    String output_vcf = mh.output_vcf
    String output_vcf_index = mh.output_vcf_index
    Int output_vcf_moved=mh.output_vcf_moved
  }
}

task MoveVcf {
  input {
    String genoset
    String GRID
    String? project_id

    String input_vcf
    String input_vcf_index

    String target_bucket
  }

  String gcs_output_dir = sub(target_bucket, "/+$", "")
  String target_folder = "~{gcs_output_dir}/~{genoset}/~{GRID}"

  String new_vcf = "${target_folder}/${basename(input_vcf)}"
  String new_vcf_index = "${target_folder}/${basename(input_vcf_index)}"

  command <<<
set +e

result=$(gsutil -q stat ~{input_vcf} || echo 1)
if [[ $result != 1 ]]; then
  echo "Source vcf file exists, moving to target bucket ..."

  set -e
    
  gsutil -m ~{"-u " + project_id} mv ~{input_vcf} \
    ~{input_vcf_index} \
    ~{target_folder}/

else
  echo "Source vcf file does not exist, checking target vcf file ..."

  result=$(gsutil -q stat ~{new_vcf} || echo 1)
  if [[ $result != 1 ]]; then
    echo "Target vcf file exists, return"
    exit 0
  else
    echo "Target vcf file does not exist, error"
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
    String output_vcf = "~{new_vcf}"
    String output_vcf_index = "~{new_vcf_index}"
    Int output_vcf_moved = 1
  }
}
