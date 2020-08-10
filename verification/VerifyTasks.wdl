version 1.0

task CompareVcfs {
  input {
    File file1
    File file2
  }

  command {
    cmp <(gunzip -c -f ~{file1} | grep -v '^#') <(gunzip -c -f ~{file2} | grep -v '^#')
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
    memory: "3 GiB"
    preemptible: 3
  }
}

task CompareGtcs {
  input {
    File file1
    File file2
    File bead_pool_manifest_file
  }

  command {
    java -Xms3g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
      CompareGtcFiles \
      --INPUT ~{file1} \
      --INPUT ~{file2} \
      --BPM_FILE ~{bead_pool_manifest_file}
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk 10 HDD"
    memory: "3.5 GiB"
    preemptible: 3
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
        diff $a $b
        if [ $? -ne 0 ]; then
            exit_code=1
            echo "Error: Files $a and $b differ" >&2
        fi
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
  }

  command {
    # get the offset of the first alignment
    test_offset="$(zcat ~{test_crai} | cut -f4 | head -n 1)"
    truth_offset="$(zcat ~{truth_crai} | cut -f4 | head -n 1)"

    # compare files with byte offset
    cmp -i "$test_offset:$truth_offset" ~{test_cram} ~{truth_cram}
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 150 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
}

task CompareCrais {

  input {
    File test_crai
    File truth_crai
  }

  command {
    # compare columns 1,2,3,5, and 6. Cannot compare column 4
    # because it is a byte offset number that may differ.
    cmp <(zcat ~{test_crai} | cut -f1,2,3,5,6) <(zcat ~{truth_crai} | cut -f1,2,3,5,6)
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
}