version 1.0

workflow VUMCMoveFile {
  input {
    String source_file
    String target_bucket
  }

  call MoveFile {
    input:
      source_file = source_file,
      target_bucket = target_bucket
  }

  output {
    String target_file = MoveFile.target_file
  }
}

task MoveFile {
  input {
    String source_file
    String target_bucket
  }

  String target_url = "${target_bucket}/${basename(source_file)}"

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

move_file ~{source_file} ~{target_url}

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_file = "~{target_url}"
  }
}