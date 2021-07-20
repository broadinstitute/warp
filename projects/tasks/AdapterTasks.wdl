version 1.0

task GetMetadata {
  input {

  }

  command {
  }

  runtime {
    docker:
    cpu: 1
    memory: "~{memory} GiB"
    disks: "local-disk ~{disk} HDD"
  }
}


