version 1.0

workflow VUMCPlink2 {
  input {
    File source_bed
    File source_bim
    File source_fam

    String plink_option

    String? parameter_file1_arg
    File? parameter_file1

    String? parameter_file2_arg
    File? parameter_file2

    String? parameter_file3_arg
    File? parameter_file3

    Array[String] suffix_list
    String target_prefix

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  scatter(suffix in suffix_list){
    String expect_file = target_prefix + suffix
  }

  call Plink2 {
    input:
      source_bed = source_bed,
      source_bim = source_bim,
      source_fam = source_fam,

      plink_option = plink_option,

      parameter_file1_arg = parameter_file1_arg,
      parameter_file1 = parameter_file1,

      parameter_file2_arg = parameter_file2_arg,
      parameter_file2 = parameter_file2,

      parameter_file3_arg = parameter_file3_arg,
      parameter_file3 = parameter_file3,

      target_prefix = target_prefix,

      expected_files = expect_file,

      docker = docker
  }

  output {
    Array[File] output_files = Plink2.output_files
  }
}

task Plink2 {
  input {
    File source_bed
    File source_bim
    File source_fam

    String plink_option

    String? parameter_file1_arg
    File? parameter_file1

    String? parameter_file2_arg
    File? parameter_file2

    String? parameter_file3_arg
    File? parameter_file3

    String target_prefix

    Array[String] expected_files

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size(source_bed, "GB") * 2) + 2

  command <<<

plink2 \
  --bed ~{source_bed} \
  --bim ~{source_bim} \
  --fam ~{source_fam} \
  ~{parameter_file1_arg + " " + parameter_file1} \
  ~{parameter_file2_arg + " " + parameter_file2} \
  ~{parameter_file3_arg + " " + parameter_file3} \
  ~{plink_option} \
  --out ~{target_prefix}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    Array[File] output_files = expected_files
  }
}