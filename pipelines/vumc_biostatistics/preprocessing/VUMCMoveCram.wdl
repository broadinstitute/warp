version 1.0

workflow VUMCMoveCram {
  input {
    String genoset
    String GRID

    String output_cram
    String output_crai
    String output_cram_md5

    String target_bucket
  }

  call MoveCram as mf {
    input:
      genoset = genoset,
      GRID = GRID,

      output_cram = output_cram,
      output_crai = output_crai,
      output_cram_md5 = output_cram_md5,

      target_bucket = target_bucket
  }

  output {
    String target_output_cram = mf.target_output_cram
    String target_output_crai = mf.target_output_crai
    String target_output_cram_md5 = mf.target_output_cram_md5

    Int target_cram_moved = mf.target_cram_moved
  }
}

task MoveCram {
  input {
    String genoset
    String GRID

    String output_cram
    String output_crai
    String output_cram_md5

    String target_bucket
  }

  String new_output_cram = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram)}"
  String new_output_crai = "${target_bucket}/${genoset}/${GRID}/${basename(output_crai)}"
  String new_output_cram_md5 = "${target_bucket}/${genoset}/${GRID}/${basename(output_cram_md5)}"

  command <<<

if [[ ! gsutil mv ~{output_cram} ~{new_output_cram} ]] ; then
  exit 1
fi

gsutil mv ~{output_crai} ~{new_output_crai}
gsutil mv ~{output_cram_md5} ~{new_output_cram_md5}
>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String target_output_cram = "~{new_output_cram}"
    String target_output_crai = "~{new_output_crai}"
    String target_output_cram_md5 = "~{new_output_cram_md5}"

    Int target_cram_moved = 1
  }
}
