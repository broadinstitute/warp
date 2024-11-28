version 1.0

task string_to_array {
  input {
    String str
    String delimiter = ","
  }
  command {
    echo ~{str} | tr '~{delimiter}' '\n'
  }
  runtime {
    cpu: 1
    docker: "ubuntu:20.04"
    preemptible: 1
    disks: "local-disk 5 HDD"
    memory: "1 GiB"
  }
  output {
    Array[String] arr = read_lines(stdout())
  }
}

task count_lines {
  input {
    File input_file
  }

  Int disk_size = ceil(size(input_file, "GB")) + 2

  command <<<
# Count the number of lines in the file
wc -l ~{input_file} > line_count.txt
  >>>

  output {
    Int num_lines = read_int("line_count.txt")
  }

  runtime {
    cpu: 1
    docker: "ubuntu:20.04"
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "1 GiB"
  }
}
