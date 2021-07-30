version 1.0

task CompareAdapterFiles {
  input {
    File test_json
    File truth_json
  }
  command <<<
  set -eo pipefail
  diff "~{test_json}" "~{truth_json}"

  if [ $? -ne 0 ];
  then
   echo "Error: ${test_json} and ${truth_json}  differ"
  fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}