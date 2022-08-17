version 1.0

workflow BuildBWARef {
  input {
      String ref_name
      File reference_fasta
      File chrom_sizes_file
  }

  call BuildBWAreference {
     input:
	    ref_name = ref_name,
            reference_fasta = reference_fasta,
            chrom_sizes_file = chrom_sizes_file
	}

	output {
	    File referenceBundle = BuildBWAreference.referenceBundle
	}
}

task BuildBWAreference {
     input {
        String ref_name ## name of the tar.bz2 files without the suffix
        File reference_fasta
        File chrom_sizes_file
     }

     command <<<
        mkdir genome
        mv ~{chrom_sizes_file} genome/chrom.sizes
        if [ ${file: -3} == ".gz" ]
        then
            gunzip -c ~{reference_fasta} > genome/genome.fa
        else
            mv ~{reference_fasta} genome/genome.fa
        fi
        bwa index genome/genome.fa
        tar cvf - genome/ > ~{reference_fasta}.tar
     >>>

     runtime {
         docker: "us.gcr.io/broad-gotc-prod/bwa:1.0.0-0.7.17-1660770463"
	 memory: "96GB"
	 disks: "local-disk 100 HDD"
	 cpu: "4"
     }

     output {
     	    File referenceBundle = "~{reference_fasta}.tar"
     }
}
