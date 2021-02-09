version 1.0

import "../verification/VerifyGermlineSingleSample.wdl" as VerifyGermlineSingleSample

workflow VerifyReprocessing {

  input {
    Array[BamPair] bam_pairs

    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_cram
    File truth_crai
    File test_cram
    File test_crai

    File truth_gvcf
    File test_gvcf
  }

  scatter(pair in bam_pairs) {
    call CompareReprocessedBams {
      input:
        test_bam = pair.test_bam,
        truth_bam = pair.truth_bam
    }
  }

  call VerifyGermlineSingleSample.VerifyGermlineSingleSample {
    input:
      truth_metrics = truth_metrics,
      test_metrics = test_metrics,
      truth_cram = truth_cram,
      truth_crai = truth_crai,
      test_cram = test_cram,
      test_crai = test_crai,
      truth_gvcf = truth_gvcf,
      test_gvcf = test_gvcf
  }
  meta {
    allowNestedInputs: true
  }
}

struct BamPair {
  File test_bam
  File truth_bam
}

task CompareReprocessedBams {

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
