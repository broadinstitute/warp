version 1.0

task CompareVcfs {
  input {
    File file1
    File file2
    Boolean fail_fast = true
  }

  command {
    fail_fast_value=~{true="true" false="" fail_fast}

    echo "Results of CompareVcfs:" >&2
    echo -e "Test:\t~{file1}" >&2
    echo -e "Truth:\t~{file2}" >&2

    cmp <(gunzip -c -f ~{file1} | grep -v '^##') <(gunzip -c -f ~{file2} | grep -v '^##') > cmp_out.txt
    exit_code=$?
    if [ "$exit_code" -eq "0" ]; then
      echo -e "Pass\tVCFs do not differ" >&2
    else
      echo -e "Fail\tVCFs differ\t`cat cmp_out.txt`" >&2
      cmp <(gunzip -c -f ~{file1} | grep -v '^##' | cut -f 1-5,7-) <(gunzip -c -f ~{file2} | grep -v '^##' | cut -f 1-5,7-) >/dev/null
      if [ $? -eq 0 ]; then
        echo -e "\tNote that differences are ONLY found in the 'quality' column" >&2
      fi
    fi

    echo $exit_code>return_code.txt
    if [[ -n $fail_fast_value ]]; then exit $exit_code; else exit 0; fi
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
    memory: "3 GiB"
    preemptible: 3
  }
  output {
    Int exit_code = read_int("return_code.txt")
    File report_file = stderr()
  }
}

task CompareGtcs {
  input {
    File file1
    File file2
    File bead_pool_manifest_file
    Boolean fail_fast = true
  }

  command {
    fail_fast_value=~{true="true" false="" fail_fast}

    java -Xms4500m -Xmx4500m -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
      CompareGtcFiles \
      --INPUT ~{file1} \
      --INPUT ~{file2} \
      --BPM_FILE ~{bead_pool_manifest_file}
    exit_code=$?
    echo $exit_code>return_code.txt
    if [[ -n $fail_fast_value ]]; then exit $exit_code; else exit 0; fi
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.6"
    disks: "local-disk 10 HDD"
    memory: "5000 MiB"
    preemptible: 3
  }
  output {
    Int exit_code = read_int("return_code.txt")
    File report_file = stderr()
  }
}

task CompareTextFiles {
  input {
    Array[File] test_text_files
    Array[File] truth_text_files
  }

  command {
    exit_code=0

    test_files_length=~{length(test_text_files)}
    truth_files_length=~{length(truth_text_files)}
    if [ $test_files_length -ne $truth_files_length ]; then
      exit_code=1
      echo "Error: Different number of input files ($test_files_length vs. $truth_files_length).  This is really not OK"
    fi

    while read -r a && read -r b <&3;
    do
      echo "Comparing File $a with $b"
      diff $a $b > diffs.txt
      if [ $? -ne 0 ];
      then
        exit_code=1
        echo "Error: Files $a and $b differ" >&2
        cat diffs.txt >&2
      fi
      # catting the diffs.txt on STDOUT as that's what's expected.
      cat diffs.txt
    done < ~{write_lines(test_text_files)} 3<~{write_lines(truth_text_files)}

    exit $exit_code
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
}


task CompareCrams {

  input {
    File test_cram
    File test_crai
    File truth_cram
    File truth_crai
    Boolean fail_fast = true
  }

  command {
    fail_fast_value=~{true="true" false="" fail_fast}

    echo "Results of CompareCrams:"
    echo -e "Test:\t~{test_cram}"
    echo -e "Truth:\t~{truth_cram}"

    # get the offset of the first alignment
    test_offset="$(zcat ~{test_crai} | cut -f4 | head -n 1)"
    truth_offset="$(zcat ~{truth_crai} | cut -f4 | head -n 1)"

    # compare files with byte offset
    cmp -i "$test_offset:$truth_offset" ~{test_cram} ~{truth_cram} > cmp_out.txt
    exit_code=$?
    if [ "$exit_code" -eq "0" ]; then
      echo -e "Pass\tCRAMs do not differ"
    else
      echo -e "Fail\tCRAMs differ\t`cat cmp_out.txt`"
    fi
    echo $exit_code>return_code.txt
    if [[ -n $fail_fast_value ]]; then exit $exit_code; else exit 0; fi
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 150 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Int exit_code = read_int("return_code.txt")
    File report_file = stdout()
  }
}

task CompareCrais {

  input {
    File test_crai
    File truth_crai
    Boolean fail_fast = true
  }

  command {
    fail_fast_value=~{true="true" false="" fail_fast}

    echo "Results of CompareCrais:"
    echo -e "Test:\t~{test_crai}"
    echo -e "Truth:\t~{truth_crai}"

    # compare columns 1,2,3,5, and 6. Cannot compare column 4
    # because it is a byte offset number that may differ.
    cmp <(zcat ~{test_crai} | cut -f1,2,3,5,6) <(zcat ~{truth_crai} | cut -f1,2,3,5,6) > cmp_out.txt
    exit_code=$?
    if [ "$exit_code" -eq "0" ]; then
      echo -e "Pass\tCRAIs do not differ"
    else
      echo -e "Fail\tCRAMs differ\t`cat cmp_out.txt`"
    fi
    echo $exit_code>return_code.txt
    if [[ -n $fail_fast_value ]]; then exit $exit_code; else exit 0; fi
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Int exit_code = read_int("return_code.txt")
    File report_file = stdout()
  }
}

task CompareGvcfs {
  input {
    File test_gvcf
    File truth_gvcf
    Boolean fail_fast = true
  }

  command {
    fail_fast_value=~{true="true" false="" fail_fast}

    echo "Results of CompareGvcfs:"
    echo -e "Test:\t~{test_gvcf}"
    echo -e "Truth:\t~{truth_gvcf}"
    exit_code=0
    DIFF_LINES=$(diff <(gunzip -c -f ~{test_gvcf} | grep -v '^##') <(gunzip -c -f ~{truth_gvcf} | grep -v '^##') | grep -e "^<" | wc -l)
    if [ $DIFF_LINES -eq 0 ]; then
      echo -e "Pass\tGVCFs do not differ"
    elif [ $DIFF_LINES -lt 10 ]; then
      echo -e "Pass\tGVCFs differ by $DIFF_LINES lines"
    else
      exit_code=1
      echo "Fail\tGVCFs differ by $DIFF_LINES lines"
      DIFF_LINES=$(diff <(gunzip -c -f ~{test_gvcf} | grep -v '^##' | cut -f 1-5,7-) <(gunzip -c -f ~{truth_gvcf} | grep -v '^##' | cut -f 1-5,7-) | grep -e "^<" | wc -l)
      if [ $DIFF_LINES -eq 0 ]; then
        echo -e "\tNote that differences are ONLY found in the 'quality' column"
      fi
    fi

    echo $exit_code>return_code.txt
    if [[ -n $fail_fast_value ]]; then exit $exit_code; else exit 0; fi
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Int exit_code = read_int("return_code.txt")
    File report_file = stdout()
  }
}

task CompareBams {

  input {
    File test_bam
    File truth_bam
    Boolean lenient_header = false
  }

  Float bam_size = size(test_bam, "GiB") + size(truth_bam, "GiB")
  Int disk_size = ceil(bam_size * 4) + 20

  command {
    set -e
    set -o pipefail

    java -Xms3500m -Xmx7000m -jar /usr/picard/picard.jar \
    CompareSAMs \
          ~{test_bam} \
          ~{truth_bam} \
          O=comparison.tsv \
          LENIENT_HEADER=~{lenient_header}
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.6"
    disks: "local-disk " + disk_size + " HDD"
    cpu: 2
    memory: "7500 MiB"
    preemptible: 3
  }
}

task CompareCompressedTextFiles {

  input {
    File test_zip
    File truth_zip
  }

  Float file_size = size(test_zip, "GiB") + size(truth_zip, "GiB")
  Int disk_size = ceil(file_size * 4) + 20

  command {
    diff <(gunzip -c -f ~{test_zip}) <(gunzip -c -f ~{truth_zip})
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
}