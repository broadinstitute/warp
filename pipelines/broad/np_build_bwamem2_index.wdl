version 1.0

workflow BuildIndex
{
  input {
    File genome_fa
    # Genome source can be NCBI or GENCODE
    String genome_source
    # Genome build is the assembly accession (NCBI) or version (GENCODE)
    String genome_build
    # Organism can be Macaque, Mouse, Human, etc.
    String organism
  }

  call BuildBWAreference {
    input:
      genome_fa = genome_fa,
      genome_source = genome_source,
      genome_build = genome_build,
      organism = organism
  }
}


task BuildBWAreference {
  input {
    File genome_fa

    # Genome source can be NCBI or GENCODE
    String genome_source
    # Genome build is the assembly accession (NCBI) or version (GENCODE)
    String genome_build
    # Organism can be Macaque, Mouse, Human, etc.
    String organism
  }

String reference_name = "bwa-mem2-2.2.1-~{organism}-~{genome_source}-build-~{genome_build}"

  command <<<
    mkdir genome
    file=~{genome_fa}
    mv ~{genome_fa} genome/genome.fa
    bwa-mem2 index genome/genome.fa
    tar --dereference -cvf - genome/ > ~{reference_name}.tar
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-bwa-mem-2:1.0.0-2.2.1_x64-linux-1685469504"
    memory: "96GB"
    disks: "local-disk 100 HDD"
    disk: "100 GB" # TES
    cpu: "4"
  }

  output {
    File reference_bundle = "~{reference_name}.tar"
  }
}
