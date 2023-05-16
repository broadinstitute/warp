version 1.0

## Copyright VUMC, 2023
##
## This WDL pipeline calculates the md5sum of a file
##
workflow VUMCMd5 {
  input {
    File input_file
  }

  call Md5File {
    input:
      input_file = input_file,
  }

  output {
    String md5sum = Md5File.md5sum
  }
}

task Md5File {
  input {
    File input_file
  }

  String md5_file = "md5.txt"
  Int disk_size = ceil(size(input_file, "GB")) + 10

  command <<<
md5sum ~{input_file} | cut -d ' ' -f1 > ~{md5_file}
>>>

  runtime {
    docker: "ubuntu:latest"
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    String md5sum = read_string("~{md5_file}")
  }
}