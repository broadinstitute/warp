version 1.0

#import "PrepareMvSusieInput.wdl" as MvSusieInput

workflow TestPrepareMvSusieInput {
  input {
    String window_id
    File window_manifest
    File window_phenotypes_tsv
    File dosage
    File phenotype_data
    File expression_covariates
    File splicing_covariates
    File protein_covariates
    File? keep_samples
    Float min_genotype_variance = 0.00000001
    Float min_phenotype_variance = 0.00000001
    String docker_image
  }

  call PrepareMvSusieInput as Prepare {
    input:
      window_id = window_id,
      window_manifest = window_manifest,
      window_phenotypes_tsv = window_phenotypes_tsv,
      dosage = dosage,
      phenotype_data = phenotype_data,
      expression_covariates = expression_covariates,
      splicing_covariates = splicing_covariates,
      protein_covariates = protein_covariates,
      keep_samples = keep_samples,
      min_genotype_variance = min_genotype_variance,
      min_phenotype_variance = min_phenotype_variance,
      docker_image = docker_image
  }

  output {
    File prepared_window = Prepare.prepared_window
    File covariate_provenance = Prepare.covariate_provenance
    File preparation_stdout = Prepare.preparation_stdout
    File preparation_stderr = Prepare.preparation_stderr
  }
}



task PrepareMvSusieInput {
  input {
    String window_id
    File window_manifest
    File window_phenotypes_tsv
    File dosage
    File phenotype_data
    File expression_covariates
    File splicing_covariates
    File protein_covariates
    File? keep_samples
    Float min_genotype_variance
    Float min_phenotype_variance
    String docker_image
  }

  command <<<
    set -euo pipefail
    output_prefix='~{window_id}'
    exec > >(tee "${output_prefix}.preparation.stdout.log") \
      2> >(tee "${output_prefix}.preparation.stderr.log" >&2)
    log() {
      printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
    }

    log "Starting preparation for ~{window_id}. The task has 16 GiB of memory."
    Rscript /opt/mvsusie/scripts/prepare_window.R \
      --windows ~{window_manifest} \
      --window-phenotypes ~{window_phenotypes_tsv} \
      --window-id ~{window_id} \
      --dosage ~{dosage} \
      --phenotype-files ~{phenotype_data} \
      --expression-covariates ~{expression_covariates} \
      --splicing-covariates ~{splicing_covariates} \
      --protein-covariates ~{protein_covariates} \
      ~{if defined(keep_samples) then "--keep-samples " + select_first([keep_samples]) else ""} \
      --min-genotype-variance ~{min_genotype_variance} \
      --min-phenotype-variance ~{min_phenotype_variance} \
      --covariate-provenance-output "${output_prefix}.preparation_covariate_provenance.tsv.gz" \
      --output "${output_prefix}.prepared_window.rds"
    log "Verifying the prepared window for ~{window_id}."
    test -s "${output_prefix}.prepared_window.rds"
    test -s "${output_prefix}.preparation_covariate_provenance.tsv.gz"
    log "Completed preparation for ~{window_id}."
  >>>

  output {
    File prepared_window = window_id + ".prepared_window.rds"
    File covariate_provenance = window_id + ".preparation_covariate_provenance.tsv.gz"
    File preparation_stdout = window_id + ".preparation.stdout.log"
    File preparation_stderr = window_id + ".preparation.stderr.log"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
    bootDiskSizeGb: 30
  }
}