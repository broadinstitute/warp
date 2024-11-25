version 1.0

task MoveOrCopyVcfFile {
  input {
    String input_vcf
    String input_vcf_index

    Boolean is_move_file = false

    String? project_id
    String target_bucket
    String genoset
    String? GRID
  }

  String action = if (is_move_file) then "mv" else "cp"

  String gcs_output_dir = sub(target_bucket, "/+$", "")

  String target_folder = if(defined(GRID)) then "~{gcs_output_dir}/~{genoset}/~{GRID}" else "~{gcs_output_dir}/~{genoset}"
  String new_vcf = "~{target_folder}/~{basename(input_vcf)}"
  String new_vcf_index = "~{target_folder}/~{basename(input_vcf_index)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} ~{input_vcf} \
  ~{input_vcf_index} \
  ~{target_folder}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String output_vcf = new_vcf
    String output_vcf_index = new_vcf_index
  }
}

task MoveOrCopyPlinkFile {
  input {
    String source_bed
    String source_bim
    String source_fam

    Boolean is_move_file = false

    String? project_id
    String target_bucket
  }

  String action = if (is_move_file) then "mv" else "cp"

  String gcs_output_dir = sub(target_bucket, "/+$", "")

  String new_bed = "~{gcs_output_dir}/~{basename(source_bed)}"
  String new_bim = "~{gcs_output_dir}/~{basename(source_bim)}"
  String new_fam = "~{gcs_output_dir}/~{basename(source_fam)}"

  command <<<

set -e

gsutil -m ~{"-u " + project_id} ~{action} ~{source_bed} \
  ~{source_bim} \
  ~{source_fam} \
  ~{gcs_output_dir}/

>>>

  runtime {
    docker: "google/cloud-sdk"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    String output_bed = new_bed
    String output_bim = new_bim
    String output_fam = new_fam
  }
}

task ExtractPgenSamples {
  input {
    File source_pgen
    File source_pvar
    File source_psam
    File extract_sample
    String chromosome

    String plink2_filter_option

    Int memory_gb = 20

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([source_pgen, source_psam, source_pvar], "GB")  * 2) + 20

  String new_pgen = chromosome + ".pgen"
  String new_pvar = chromosome + ".pvar"
  String new_psam = chromosome + ".psam"

  command <<<

plink2 \
  --pgen ~{source_pgen} \
  --pvar ~{source_pvar} \
  --psam ~{source_psam} \
  ~{plink2_filter_option} \
  --keep ~{extract_sample} \
  --make-pgen \
  --out ~{chromosome}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_pgen = new_pgen
    File output_pvar = new_pvar
    File output_psam = new_psam
  }
}


task ExtractPgenRegions {
  input {
    File source_pgen
    File source_pvar
    File source_psam
    File region_bed
    String chromosome

    String plink2_filter_option

    Int memory_gb = 20

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([source_pgen, source_psam, source_pvar], "GB")  * 2) + 20

  String new_pgen = chromosome + ".pgen"
  String new_pvar = chromosome + ".pvar"
  String new_psam = chromosome + ".psam"

  command <<<

plink2 ~{plink2_filter_option} \
  --pgen ~{source_pgen} \
  --pvar ~{source_pvar} \
  --psam ~{source_psam} \
  --extract bed0 ~{region_bed} \
  --make-pgen \
  --out ~{chromosome}

>>>

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_pgen = new_pgen
    File output_pvar = new_pvar
    File output_psam = new_psam
  }
}

task MergePgenFiles {
  input {
    Array[File] pgen_files
    Array[File] pvar_files
    Array[File] psam_files

    String output_prefix

    Int memory_gb = 20
    Int cpu = 8

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil((size(pgen_files, "GB") + size(pvar_files, "GB") + size(psam_files, "GB"))  * 3) + 20

  String new_pgen = output_prefix + ".pgen"
  String new_pvar = output_prefix + ".pvar"
  String new_psam = output_prefix + ".psam"

  String new_merged_pgen = output_prefix + "-merge.pgen"
  String new_merged_pvar = output_prefix + "-merge.pvar"
  String new_merged_psam = output_prefix + "-merge.psam"

  command <<<

cat ~{write_lines(pgen_files)} > pgen.list
cat ~{write_lines(pvar_files)} > pvar.list
cat ~{write_lines(psam_files)} > psam.list

paste pgen.list pvar.list psam.list > merge.list

plink2 --pmerge-list merge.list --make-pgen --out ~{output_prefix} --threads ~{cpu}

rm -f ~{new_pgen} ~{new_pvar} ~{new_psam}

mv ~{new_merged_pgen} ~{new_pgen}
mv ~{new_merged_pvar} ~{new_pvar}
mv ~{new_merged_psam} ~{new_psam}

>>>

  runtime {
    cpu: cpu
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_pgen = new_pgen
    File output_pvar = new_pvar
    File output_psam = new_psam
  }
}
