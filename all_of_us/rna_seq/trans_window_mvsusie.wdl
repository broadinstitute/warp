version 1.0

workflow TransWindowMvSusie {
  input {
    String window_id
    File? prepared_window
    File? window_manifest
    File? window_phenotypes_tsv
    File? dosage
    File? phenotype_data
    File? expression_covariates
    File? splicing_covariates
    File? protein_covariates
    File? keep_samples
    Int start_L = 10
    Int step_L = 5
    Int max_L = 40
    Float greedy_lbf_cutoff = 1.0
    Int max_iter = 100
    Float tol = 0.0001
    Float coverage = 0.95
    Float min_abs_corr = 0.5
    Float min_genotype_variance = 0.00000001
    Float min_phenotype_variance = 0.00000001
    Int n_thread = 1
    Int mashr_n_pca = 5
    Float mashr_strong_lfsr = 0.05
    Int? mashr_seed
    String docker_image = "ghcr.io/aou-multiomics-analysis/mvsusier-trans-window-mvsusie:latest"
  }

  call ValidateMvSusieInputs {
    input:
      window_id = window_id,
      has_prepared_window = defined(prepared_window),
      has_window_manifest = defined(window_manifest),
      has_window_phenotypes_tsv = defined(window_phenotypes_tsv),
      has_dosage = defined(dosage),
      has_phenotype_data = defined(phenotype_data),
      has_expression_covariates = defined(expression_covariates),
      has_splicing_covariates = defined(splicing_covariates),
      has_protein_covariates = defined(protein_covariates),
      docker_image = docker_image
  }

  if (ValidateMvSusieInputs.run_preparation) {
    call PrepareMvSusieInput {
      input:
        window_id = ValidateMvSusieInputs.validated_window_id,
        window_manifest = select_first([window_manifest]),
        window_phenotypes_tsv = select_first([window_phenotypes_tsv]),
        dosage = select_first([dosage]),
        phenotype_data = select_first([phenotype_data]),
        expression_covariates = select_first([expression_covariates]),
        splicing_covariates = select_first([splicing_covariates]),
        protein_covariates = select_first([protein_covariates]),
        keep_samples = keep_samples,
        min_genotype_variance = min_genotype_variance,
        min_phenotype_variance = min_phenotype_variance,
        docker_image = docker_image
    }
  }

  File resolved_prepared_window = select_first([prepared_window, PrepareMvSusieInput.prepared_window])

  call FitMvSusie {
    input:
      window_id = ValidateMvSusieInputs.validated_window_id,
      prepared_window = resolved_prepared_window,
      start_L = start_L,
      step_L = step_L,
      max_L = max_L,
      greedy_lbf_cutoff = greedy_lbf_cutoff,
      max_iter = max_iter,
      tol = tol,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      n_thread = n_thread,
      mashr_n_pca = mashr_n_pca,
      mashr_strong_lfsr = mashr_strong_lfsr,
      mashr_seed = mashr_seed,
      docker_image = docker_image
  }

  call SummarizeMvSusie {
    input:
      window_id = ValidateMvSusieInputs.validated_window_id,
      prepared_window = FitMvSusie.prepared_window_output,
      mvsusie_fit = FitMvSusie.mvsusie_fit,
      docker_image = docker_image
  }

  call PlotMvSusie {
    input:
      window_id = ValidateMvSusieInputs.validated_window_id,
      prepared_window = FitMvSusie.prepared_window_output,
      mvsusie_fit = FitMvSusie.mvsusie_fit,
      docker_image = docker_image
  }

  output {
    File prepared_window_output = FitMvSusie.prepared_window_output
    File mvsusie_fit = FitMvSusie.mvsusie_fit
    File mashr_training = FitMvSusie.mashr_training
    File greedy_L_history = FitMvSusie.greedy_L_history
    File covariate_provenance = FitMvSusie.covariate_provenance
    File run_stdout = FitMvSusie.run_stdout
    File run_stderr = FitMvSusie.run_stderr
    File session_info = FitMvSusie.session_info
    File variant_pip = SummarizeMvSusie.variant_pip
    File credible_sets = SummarizeMvSusie.credible_sets
    File credible_set_members = SummarizeMvSusie.credible_set_members
    File component_feature_support = SummarizeMvSusie.component_feature_support
    File window_qc = SummarizeMvSusie.window_qc
    File effect_plot_png = PlotMvSusie.effect_plot_png
    File effect_plot_pdf = PlotMvSusie.effect_plot_pdf
    File effect_plot_rds = PlotMvSusie.effect_plot_rds
  }
}

task ValidateMvSusieInputs {
  input {
    String window_id
    Boolean has_prepared_window
    Boolean has_window_manifest
    Boolean has_window_phenotypes_tsv
    Boolean has_dosage
    Boolean has_phenotype_data
    Boolean has_expression_covariates
    Boolean has_splicing_covariates
    Boolean has_protein_covariates
    String docker_image
  }

  command <<<
    set -euo pipefail
    log() {
      printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
    }

    window_id='~{window_id}'
    if [[ ! "$window_id" =~ ^[A-Za-z0-9._-]+$ ]]; then
      log "window_id may contain letters, numbers, periods, underscores, and hyphens only."
      exit 1
    fi
    printf '%s\n' "$window_id" > validated_window_id.txt

    if ~{has_prepared_window}; then
      log "A prepared window is present. The workflow will skip phenotype and genotype preparation."
      printf 'false\n' > run_preparation.txt
      exit 0
    fi

    log "No prepared window is present. Checking the raw preparation inputs."
    missing=()
    if ! ~{has_window_manifest}; then missing+=("window_manifest"); fi
    if ! ~{has_window_phenotypes_tsv}; then missing+=("window_phenotypes_tsv"); fi
    if ! ~{has_dosage}; then missing+=("dosage"); fi
    if ! ~{has_phenotype_data}; then missing+=("phenotype_data"); fi
    if ! ~{has_expression_covariates}; then missing+=("expression_covariates"); fi
    if ! ~{has_splicing_covariates}; then missing+=("splicing_covariates"); fi
    if ! ~{has_protein_covariates}; then missing+=("protein_covariates"); fi

    if (( ${#missing[@]} > 0 )); then
      log "Raw preparation inputs are missing: ${missing[*]}."
      exit 1
    fi

    log "All raw preparation inputs are present."
    printf 'true\n' > run_preparation.txt
  >>>

  output {
    Boolean run_preparation = read_boolean("run_preparation.txt")
    String validated_window_id = read_string("validated_window_id.txt")
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "1 GiB"
    disks: "local-disk 10 SSD"
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
  }
}

task FitMvSusie {
  input {
    String window_id
    File prepared_window
    Int start_L
    Int step_L
    Int max_L
    Float greedy_lbf_cutoff
    Int max_iter
    Float tol
    Float coverage
    Float min_abs_corr
    Int n_thread
    Int mashr_n_pca
    Float mashr_strong_lfsr
    Int? mashr_seed
    String docker_image
  }

  command <<<
    set -euo pipefail
    output_prefix='~{window_id}'
    exec > >(tee "${output_prefix}.run.stdout.log") \
      2> >(tee "${output_prefix}.run.stderr.log" >&2)
    log() {
      printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
    }

    cp '~{prepared_window}' "${output_prefix}.prepared_window.rds"
    log "Starting joint mvSuSiE for ~{window_id}. The task has 8 GiB of memory."
    log "Greedy L starts at ~{start_L}, increases by ~{step_L}, and stops at ~{max_L}; the cutoff is ~{greedy_lbf_cutoff}."
    log "The fit is verbose. Iteration updates will be present in the task log."
    Rscript /opt/mvsusie/scripts/fit_window.R \
      --prepared "${output_prefix}.prepared_window.rds" \
      --window-id ~{window_id} \
      --start-L ~{start_L} \
      --step-L ~{step_L} \
      --max-L ~{max_L} \
      --greedy-lbf-cutoff ~{greedy_lbf_cutoff} \
      --max-iter ~{max_iter} \
      --tol ~{tol} \
      --coverage ~{coverage} \
      --min-abs-corr ~{min_abs_corr} \
      --n-thread ~{n_thread} \
      --mashr-n-pca ~{mashr_n_pca} \
      --mashr-strong-lfsr ~{mashr_strong_lfsr} \
      ~{if defined(mashr_seed) then "--mashr-seed " + select_first([mashr_seed]) else ""} \
      --covariate-provenance-output "${output_prefix}.covariate_provenance.tsv.gz" \
      --mashr-output "${output_prefix}.mashr_training.rds" \
      --greedy-history-output "${output_prefix}.greedy_L_history.tsv" \
      --output "${output_prefix}.mvsusie_fit.rds"
    log "Writing the R session information."
    Rscript -e 'writeLines(capture.output(sessionInfo()), commandArgs(TRUE)[[1L]])' \
      "${output_prefix}.session_info.txt"
    log "Verifying the joint model outputs for ~{window_id}."
    test -s "${output_prefix}.prepared_window.rds"
    test -s "${output_prefix}.mvsusie_fit.rds"
    test -s "${output_prefix}.mashr_training.rds"
    test -s "${output_prefix}.greedy_L_history.tsv"
    test -s "${output_prefix}.covariate_provenance.tsv.gz"
    test -s "${output_prefix}.session_info.txt"
    log "Completed joint mvSuSiE for ~{window_id}."
  >>>

  output {
    File prepared_window_output = window_id + ".prepared_window.rds"
    File mvsusie_fit = window_id + ".mvsusie_fit.rds"
    File mashr_training = window_id + ".mashr_training.rds"
    File greedy_L_history = window_id + ".greedy_L_history.tsv"
    File covariate_provenance = window_id + ".covariate_provenance.tsv.gz"
    File run_stdout = window_id + ".run.stdout.log"
    File run_stderr = window_id + ".run.stderr.log"
    File session_info = window_id + ".session_info.txt"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "8 GiB"
    disks: "local-disk 500 SSD"
  }
}

task SummarizeMvSusie {
  input {
    String window_id
    File prepared_window
    File mvsusie_fit
    String docker_image
  }

  command <<<
    set -euo pipefail
    output_prefix='~{window_id}'
    log() {
      printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
    }
    log "Starting the joint mvSuSiE summary."
    mkdir -p window_outputs
    Rscript /opt/mvsusie/scripts/summarize_window.R \
      --prepared ~{prepared_window} \
      --fit ~{mvsusie_fit} \
      --output-dir window_outputs
    mv window_outputs/variant_pip.tsv.gz "${output_prefix}.variant_pip.tsv.gz"
    mv window_outputs/credible_sets.tsv.gz "${output_prefix}.credible_sets.tsv.gz"
    mv window_outputs/credible_set_members.tsv.gz "${output_prefix}.credible_set_members.tsv.gz"
    mv window_outputs/component_feature_support.tsv.gz "${output_prefix}.component_feature_support.tsv.gz"
    mv window_outputs/window_qc.tsv "${output_prefix}.window_qc.tsv"
    log "Verifying the joint mvSuSiE summary outputs."
    test -s "${output_prefix}.variant_pip.tsv.gz"
    test -e "${output_prefix}.credible_sets.tsv.gz"
    test -e "${output_prefix}.credible_set_members.tsv.gz"
    test -s "${output_prefix}.component_feature_support.tsv.gz"
    test -s "${output_prefix}.window_qc.tsv"
    log "Completed the joint mvSuSiE summary."
  >>>

  output {
    File variant_pip = window_id + ".variant_pip.tsv.gz"
    File credible_sets = window_id + ".credible_sets.tsv.gz"
    File credible_set_members = window_id + ".credible_set_members.tsv.gz"
    File component_feature_support = window_id + ".component_feature_support.tsv.gz"
    File window_qc = window_id + ".window_qc.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}

task PlotMvSusie {
  input {
    String window_id
    File prepared_window
    File mvsusie_fit
    String docker_image
  }

  command <<<
    set -euo pipefail
    output_prefix='~{window_id}'
    log() {
      printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
    }
    log "Starting the mvSuSiE credible-set-by-feature plot."
    Rscript /opt/mvsusie/scripts/plot_window_mvsusie.R \
      --prepared ~{prepared_window} \
      --fit ~{mvsusie_fit} \
      --png "${output_prefix}.effect_plot.png" \
      --pdf "${output_prefix}.effect_plot.pdf" \
      --plot-rds "${output_prefix}.effect_plot.rds"
    log "Verifying the mvSuSiE plot outputs."
    test -s "${output_prefix}.effect_plot.png"
    test -s "${output_prefix}.effect_plot.pdf"
    test -s "${output_prefix}.effect_plot.rds"
    log "Completed the mvSuSiE credible-set-by-feature plot."
  >>>

  output {
    File effect_plot_png = window_id + ".effect_plot.png"
    File effect_plot_pdf = window_id + ".effect_plot.pdf"
    File effect_plot_rds = window_id + ".effect_plot.rds"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}
