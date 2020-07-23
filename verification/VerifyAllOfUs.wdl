version 1.0

import "../verification/VerifyGermlineSingleSample.wdl"

workflow VerifyAllOfUs {
  input {
    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_cram
    File truth_crai
    File test_cram
    File test_crai

    File truth_gvcf
    File test_gvcf

    File truth_variants
    File test_significant_variants_vcf

    File truth_filtration_report
    File test_filtration_report
  }

  call VerifyGermlineSingleSample.VerifyGermlineSingleSample as VerifyGermlineSingleSampleWorkflow {
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

  call CompareFiltrationReport {
    input:
      test_filtration_report = test_filtration_report,
      truth_filtration_report = truth_filtration_report
  }

  call VerifySignificantVariantsVcf {
    input:
      test_variants = test_significant_variants_vcf,
      truth_variants = truth_variants
  }
}

task CompareFiltrationReport {
  input {
    File test_filtration_report
    File truth_filtration_report
  }

  command {
    # The report lists vcfs in the cromwell execution bucket, so we have to parse it line by line and munge it. Argh.
    touch munged.tsv
    while read line; do
      # Take everything after the last slash and print it, line by line, to munged.tsv
      echo "$line" | sed 's:.*/::' >> munged.tsv
    done < ~{test_filtration_report}

    diff <(sort ~{truth_filtration_report}) <(sort munged.tsv) > diff.txt
    # Diff should be empty. Output it in case it isn't.
    test ! -s diff.txt
  }

  runtime {
    docker: "phusion/baseimage"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }

  output {
    File diff = "diff.txt"
  }
}

task VerifySignificantVariantsVcf {

  input {
    File test_variants
    File truth_variants
  }

  command <<<
    gunzip -c ~{test_variants} | grep -v '^#' | awk '{print $1","$2","$4","$5}' > found_variants.csv

    while IFS=$',' read chromosome position rsid ref alt reason ; do
        if grep "${chromosome},${position},${ref},${alt}" found_variants.csv > /dev/null; then
            echo "Found mutation ${chromosome}:${position} ${ref}->${alt}"
        else
            echo "Missing mutation ${chromosome}:${position} ${ref}->${alt}"
            RETURN=1
        fi
    done < ~{truth_variants}

    exit ${RETURN}
  >>>

  runtime {
    docker: "phusion/baseimage"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }

  output {
    File found_variants = "found_variants.csv"
  }
}
