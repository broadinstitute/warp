version 1.0

# testing ToA and write_lines([]) as task input
workflow test_empty_write_lines_input {

  input {
    File write_lines_at_wdl_input = write_lines([])
  }

  # use file defined at wdl input
  call LocalizeFile as LocalizeEmptyFileFromWdlInput {
    input:
      input_file = write_lines_at_wdl_input
  }

  # use file defined at task input
  call LocalizeFile as LocalizeEmptyFileFromTaskInput {
    input:
      input_file = write_lines([])
  }
}

task LocalizeFile {
  input {
    File input_file
  }
    
  command <<<
    set -euo pipefail

    cat ~{input_file} | wc -l > num_lines.txt
  >>>

  runtime {
    docker: "ubuntu:20.04"
    memory: "2 GiB"
    cpu: "1"
    disks: "local-disk 8 HDD"
  }

  output {
    Int num_lines = read_int("num_lines.txt")
  }
}
