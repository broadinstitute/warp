version 1.0

task ScatterBam {
  input {
    File bam_to_scatter
    Int scatter_width
  }
  
  Int disk_size = ceil(size(bam_to_scatter, "GiB") * 3)

  command {
    mkdir scattered_bams
    java -Xms7g -jar /usr/picard/picard.jar \
      SplitSamByNumberOfReads \
      INPUT=${bam_to_scatter} \
      SPLIT_TO_N_FILES=${scatter_width} \
      OUT_PREFIX=${basename(bam_to_scatter, '.bam')}_split \
      OUTPUT=scattered_bams
  }

  output {
    Array[File] scattered_bams = glob("scattered_bams/*")
  }

  runtime {
    disks: "local-disk ${disk_size} HDD"
    cpu: 2
    memory: "7.5 GiB"
    docker: "quay.io/humancellatlas/secondary-analysis-picard:2.20.4"
  }
}
