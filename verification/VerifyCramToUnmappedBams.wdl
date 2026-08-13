version 1.0


workflow VerifyCramToUnmappedBams {

  input {
    Array[File] truth_bam
    Array[File] test_bam
    Boolean? done
  }

  scatter(idx in range(length(truth_bam))) {
    call CompareBams {
      input:
        test_bam = test_bam[idx],
        truth_bam = truth_bam[idx]
    }
  }

  meta {
    allowNestedInputs: true
  }
}

task CompareBams {

  input {
    File test_bam
    File truth_bam
  }

  Float bam_size = size(test_bam, "GiB") + size(truth_bam, "GiB")
  Int disk_size = ceil(bam_size * 4) + 20

  command {
    # ApplyBQSR changes quality scores, so we need to strip those out of both BAMs
    cmp \
      <(samtools sort -n ~{test_bam} | samtools view | cut -d$'\t' -f 1-10,12-) \
      <(samtools sort -n ~{truth_bam} | samtools view | cut -d$'\t' -f 1-10,12-)
  }

  runtime {
    docker: "biocontainers/samtools:1.3.1"
    disks: "local-disk " + disk_size + " HDD"
    cpu: 2
    memory: "7.5 GiB"
    preemptible: 3
  }

}
