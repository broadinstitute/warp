version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils

workflow VUMCGenotypePCA {
  input {
    File pgen_file
    File pvar_file
    File psam_file

    String target_prefix

    Int n_geno_pcs = 10

    String? project_id
    String? target_gcp_folder
  }

  call PlinkPCA {
    input:
      pgen_file = pgen_file,
      pvar_file = pvar_file,
      psam_file = psam_file,
      n_geno_pcs = n_geno_pcs,
      target_prefix = target_prefix
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyOneFile as CopyFile {
      input:
        source_file = PlinkPCA.output_pca_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_pca_file = select_first([CopyFile.output_file, PlinkPCA.output_pca_file])
  }
}

task PlinkPCA {
  input {
    File pgen_file
    File pvar_file
    File psam_file

    Int n_geno_pcs

    String target_prefix

    Int memory_gb = 20
    Int cpu = 8

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([pgen_file, pvar_file, psam_file], "GB")  * 2) + 20

  String pca_file = target_prefix + ".genotype.pca.txt"

  command <<<

plink2 \
  --pgen ~{pgen_file} \
  --pvar ~{pvar_file} \
  --psam ~{psam_file} \
  --indep-pairwise 50000 200 0.05 \
  --out ~{target_prefix}.pruned_variants \
  --const-fid --set-all-var-ids @:#:\$r:\$a \
  --new-id-max-allele-len 1000

plink2 \
  --pgen ~{pgen_file} \
  --pvar ~{pvar_file} \
  --psam ~{psam_file} \
  --extract ~{target_prefix}.pruned_variants.prune.in \
  --make-pgen \
  --out ~{target_prefix}.pruned \
  --const-fid \
  --set-all-var-ids @:#:\$r:\$a \
  --new-id-max-allele-len 1000

plink2 \
  -pfile ~{target_prefix}.pruned \
  --pca ~{n_geno_pcs} \
  --out ~{target_prefix}.pruned

mv ~{target_prefix}.pruned.eigenvec ~{pca_file}

>>>

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_pca_file = "~{pca_file}"
  }
}
