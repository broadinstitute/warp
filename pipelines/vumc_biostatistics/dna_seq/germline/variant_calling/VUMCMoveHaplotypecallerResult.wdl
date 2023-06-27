version 1.0

workflow VUMCMoveHaplotypecallerResult {
  input {
    String genoset
    String GRID

    String output_vcf
    String output_vcf_index

    String target_bucket
  }

  call MoveHaplotypecallerResult as mh {
    input:
      genoset = genoset,
      GRID = GRID,

      output_vcf = output_vcf,
      output_vcf_index = output_vcf_index,

      target_bucket = target_bucket
  }

  output {
    String target_output_vcf = mh.target_output_vcf
    String target_output_vcf_index = mh.target_output_vcf_index
    Int target_vcf_moved=mh.target_vcf_moved
  }
}

task MoveHaplotypecallerResult {
  input {
    String genoset
    String GRID

    String output_vcf
    String output_vcf_index

    String target_bucket
  }

  String new_output_vcf = "${target_bucket}/${genoset}/${GRID}/${basename(output_vcf)}"
  String new_output_vcf_index = "${target_bucket}/${genoset}/${GRID}/${basename(output_vcf_index)}"

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

move_file ~{output_vcf} ~{new_output_vcf}

move_file ~{output_vcf_index} ~{new_output_vcf_index}

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_output_vcf = "~{new_output_vcf}"
    String target_output_vcf_index = "~{new_output_vcf_index}"
    Int target_vcf_moved = 1
  }
}
