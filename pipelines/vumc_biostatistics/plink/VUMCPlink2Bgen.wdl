version 1.0

workflow VUMCPlink2Bgen {
  input {
    File source_bgen
    File source_sample

    String bgen_ref_alt_mode = "ref-last"

    String target_prefix

    String plink_option

    String? parameter_file1_arg
    File? parameter_file1

    String? parameter_file2_arg
    File? parameter_file2

    String? parameter_file3_arg
    File? parameter_file3

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

    Int? memory_size=10
  }

  call Plink2Bgen {
    input:
      source_bgen = source_bgen,
      source_sample = source_sample,

      bgen_ref_alt_mode = bgen_ref_alt_mode,

      target_prefix = target_prefix,

      plink_option = plink_option,

      parameter_file1_arg = parameter_file1_arg,
      parameter_file1 = parameter_file1,

      parameter_file2_arg = parameter_file2_arg,
      parameter_file2 = parameter_file2,

      parameter_file3_arg = parameter_file3_arg,
      parameter_file3 = parameter_file3,

      docker = docker,

      memory_size = memory_size
  }

  output {
    File output_bgen = Plink2Bgen.output_bgen
    File output_sample = Plink2Bgen.output_sample
  }
}

task Plink2Bgen {
  input {
    File source_bgen
    File source_sample

    String bgen_ref_alt_mode    

    String target_prefix

    String plink_option

    String? parameter_file1_arg
    File? parameter_file1

    String? parameter_file2_arg
    File? parameter_file2

    String? parameter_file3_arg
    File? parameter_file3

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

    Int? memory_size=10
  }

  Int disk_size = ceil(size(source_bgen, "GB") * 2) + 2

  String new_bgen = "${target_prefix}.bgen"
  String new_sample = "${target_prefix}.sample"

  command <<<

plink2 \
  --bgen  ~{source_bgen} ~{bgen_ref_alt_mode} \
  --sample ~{source_sample} \
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
    memory: memory_size + " GiB"
  }
  output {
    File output_bgen = new_bgen
    File output_sample = new_sample
  }
}