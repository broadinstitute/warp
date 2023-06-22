version 1.0

workflow VUMCMoveCram {
  input {
    String genoset
    String GRID

    String output_cram
    String output_cram_index
    String output_cram_md5

    String target_bucket
  }

  call MoveCram as mf {
    input:
      genoset = genoset,
      GRID = GRID,

      output_cram = output_cram,
      output_cram_index = output_cram_index,
      output_cram_md5 = output_cram_md5,

      target_bucket = target_bucket
  }

  output {
    String target_output_cram = mf.target_output_cram
    String target_output_cram_index = mf.target_output_cram_index
    String target_output_cram_md5 = mf.target_output_cram_md5

    Int target_cram_moved = mf.target_cram_moved
  }
}

task MoveCram {
  input {
    String genoset
    String GRID

    String output_cram
    String output_cram_index
    String output_cram_md5

    String target_bucket
  }

  String new_output_cram = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram)}"
  String new_output_cram_index = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram_index)}"
  String new_output_cram_md5 = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram_md5)}"

  command <<<

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

    gsutil mv $SOURCE_FILE $TARGET_FILE
    status=$?
    return $status
  fi
}

move_file ~{output_cram_index} ~{new_output_cram_index}
move_file ~{output_cram_md5} ~{new_output_cram_md5}
move_file ~{output_cram} ~{new_output_cram}
>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_output_cram = "~{new_output_cram}"
    String target_output_cram_index = "~{new_output_cram_index}"
    String target_output_cram_md5 = "~{new_output_cram_md5}"

    Int target_cram_moved = 1
  }
}
