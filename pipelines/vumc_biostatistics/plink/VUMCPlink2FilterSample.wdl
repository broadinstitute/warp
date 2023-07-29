version 1.0

workflow VUMCPlink2FilterSample {
  input {
    File source_bed
    File source_bim
    File source_fam

    File sample_file

    String target_prefix

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  call Plink2Filter {
    input:
      source_bed = source_bed,
      source_bim = source_bim,
      source_fam = source_fam,

      sample_file = sample_file,

      target_prefix = target_prefix,

      docker = docker
  }

  output {
    File target_bed = Plink2Filter.target_bed
    File target_bim = Plink2Filter.target_bim
    File target_fam = Plink2Filter.target_fam
  }
}

task Plink2Filter {
  input {
      File source_bed
      File source_bim
      File source_fam

      File sample_file

      String target_prefix

      String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size(source_bed, "GB") * 2) + 2

  String output_bed = "${target_prefix}.bed"
  String output_bim = "${target_prefix}.bim"
  String output_fam = "${target_prefix}.fam"

  command <<<

plink2 \
  --bed ~{source_bed} \
  --bim ~{source_bim} \
  --fam ~{source_fam} \
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
    File target_bed = output_bed
    File target_bim = output_bim
    File target_fam = output_fam
  }
}