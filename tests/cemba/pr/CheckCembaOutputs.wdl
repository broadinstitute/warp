workflow CheckCembaOutputs {
  File cemba_bam_output_file_read1
  File cemba_bam_output_file_read2
  File cemba_truth_file_read1
  File cemba_truth_file_read2

  call CompareBams as CompareBamsRead1 {
    input:
      cemba_bam_output_file = cemba_bam_output_file_read1,
      cemba_truth_file = cemba_truth_file_read1
  }

  call CompareBams as CompareBamsRead2 {
    input:
      cemba_bam_output_file = cemba_bam_output_file_read2,
      cemba_truth_file = cemba_truth_file_read2
  }
}

task CompareBams {
  File cemba_bam_output_file
  File cemba_truth_file

  # input file size
  Float input_size = size(cemba_bam_output_file, "GB") + size(cemba_truth_file, "GB")

  command <<<
    # compare the 2 files (stream outputs of samtools view into cmp)
    # cmp instead of diff for truncation/fail issue using diff w/ larger files
    cmp <(samtools view ${cemba_bam_output_file}) <(samtools view ${cemba_truth_file})
  >>>

  runtime {
    docker: "us.gcr.io/broad-biccn-dev/samtools:0.11"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 1 * input file size
    disks: "local-disk " + ceil(1 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }
}
