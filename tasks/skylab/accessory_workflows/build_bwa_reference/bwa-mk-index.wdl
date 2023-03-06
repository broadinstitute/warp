version 1.0

workflow BuildBWARef {
  input {
      File genome_fa
      File chrom_sizes_file
      # Organism can be Macaque, Mouse, Human, etc.
      String organism
      # Genome source can be NCBI or GENCODE
      String genome_source
  }


  call BuildBWAreference {
     input:
        genome_fa = genome_fa,
        chrom_sizes_file = chrom_sizes_file,
        organism = organism,
        genome_source = genome_source
	}

	output {
	  File reference_bundle = BuildBWAreference.reference_bundle
    }
  }

task BuildBWAreference {
     input {
        File genome_fa
        File chrom_sizes_file

        # Organism can be Macaque, Mouse, Human, etc.
        String organism
        # Genome source can be NCBI or GENCODE
        String genome_source
    }

    String reference_name = "bwa0.7.17-~{organism}-~{genome_source}"
    String reference_bundle = "~{reference_name}.tar"

     command <<<
        mkdir genome
        mv ~{chrom_sizes_file} genome/chrom.sizes
        file=~{genome_fa}
        if [ ${file: -3} == ".gz" ]
        then
            gunzip -c ~{genome_fa} > genome/genome.fa
        else
            mv ~{genome_fa} genome/genome.fa
        fi
        bwa index genome/genome.fa
        tar --dereference -cvf - genome/ > ~{reference_bundle}.tar
     >>>

     runtime {
       docker: "us.gcr.io/broad-gotc-prod/bwa:1.0.0-0.7.17-1660770463"
	   memory: "96GB"
	   disks: "local-disk 100 HDD"
       disk: "100 GB" # TES
	   cpu: "4"
     }

     output {
       File reference_bundle = "~{reference_bundle}.tar"
     }
}
