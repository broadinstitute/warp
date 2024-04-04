version 1.0

workflow TestCountVariants {

  input {
    File vcf_path
    File vcf_index_path
  }

  # Lift over the array to hg38.
  call CountVariantsTest {
    input:
      vcf = vcf_path,
      vcf_index = vcf_index_path
  }

  output {
    File test1 = CountVariantsTest.test1
    File test2 = CountVariantsTest.test2
    File test3 = CountVariantsTest.test3
  }
}

task CountVariantsTest {
  input {
    File vcf
    File vcf_index


    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int cpu = 1
    Int memory_mb = 8000
    Int disk_size_gb = 2 * ceil(size([vcf, vcf_index], "GiB")) + 20
  }
  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    set -e -o pipefail

    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" CountVariants -V ~{vcf} > test_0.txt
    cat test_0.txt | tail -n 1 > test_1.txt
    echo "test_1:"
    cat test_1.txt
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" CountVariants -V ~{vcf} | tail -n 1 > test_2.txt
    echo "test_2:"
    cat test_2.txt
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" CountVariants -V ~{vcf} 2>&1 | tail -n 1 > test_3.txt
    echo "test_3:"
    cat test_3.txt
  >>>
  output {
    File test1 = "test_1.txt"
    File test2 = "test_2.txt"
    File test3 = "test_3.txt"
  }
  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
  }
}