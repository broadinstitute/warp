version 1.0

workflow VUMCPgenFilterSample {
  input {
    File source_pgen
    File source_pvar
    File source_psam

    File sample_file

    String target_prefix

    String? plink2_option = "--chr-set 22 no-xy"

    String? project_id
    String? target_gcp_folder
  }

  call Plink2Filter {
    input:
      source_pgen = source_pgen,
      source_pvar = source_pvar,
      source_psam = source_psam,

      plink2_option = plink2_option,

      sample_file = sample_file,

      target_prefix = target_prefix,
  }

  output {
    File target_pgen = Plink2Filter.target_pgen
    File target_pvar = Plink2Filter.target_pvar
    File target_psam = Plink2Filter.target_psam
  }
}

task Plink2Filter {
  input {
      File source_pgen
      File source_pvar
      File source_psam

      String? plink2_option

      File sample_file

      String target_prefix

      String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size(source_pgen, "GB") * 2) + 2

  String output_pgen = "${target_prefix}.bed"
  String output_pvar = "${target_prefix}.bim"
  String output_psam = "${target_prefix}.fam"

  command <<<

plink2 ~{plink2_option} \
  --pgen ~{source_pgen} \
  --pvar ~{source_pvar} \
  --psam ~{source_psam} \
  --keep ~{sample_file} \
  --make-bed \
  --out ~{target_prefix}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    File target_pgen = output_pgen
    File target_pvar = output_pvar
    File target_psam = output_psam
  }
}