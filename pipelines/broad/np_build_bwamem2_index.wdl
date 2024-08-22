version 1.0

workflow BuildSeqDict
{
  input {
    File genome_fa
  }

  call picardCreateDict {
    input:
      genome_fa = genome_fa
}
}


task picardCreateDict {
  input {
    File genome_fa
  }

  command <<<
    java -jar /usr/local/bin/picard.jar CreateSequenceDictionary R=~{genome_fa}  O=Homo_sapiens_assembly38.dict
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-bwa-mem2:latest"
    memory: "96GB"
    disks: "local-disk 100 HDD"
    disk: "100 GB" # TES
    cpu: "4"
  }

  output {
    File seq_dict = "Homo_sapiens_assembly38.dict"
  }
}
