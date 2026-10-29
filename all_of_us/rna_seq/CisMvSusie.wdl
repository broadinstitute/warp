version 1.0

# Cis multivariate fine-mapping (mvSuSiE) across expression, splicing and protein.
#
# One window per gene (TSS +/- cis_window_bp). The outcomes in each window are all
# features of that gene in each modality. The collaborator's trans scripts are reused
# unchanged: BuildCisWindows writes per-gene "associations" and "targets" tables in the
# format those scripts expect, with the top-N caps set high enough to keep every feature.

workflow CisMvSusie {
  input {
    # Genotypes: bgzipped wide dosage TSV (CHROM POS REF ALT <samples>) plus .tbi
    File genome_dosage
    File genome_dosage_tbi

    # tensorQTL phenotype BEDs (#chr start end phenotype_id <samples>)
    File expression_phenotypes
    File splicing_phenotypes
    File protein_phenotypes

    # tensorQTL cis permutation results (tsv or tsv.gz), one per modality
    File expression_cis_results
    File splicing_cis_results
    File protein_cis_results

    # One covariate file per modality (covariates as rows, samples as columns)
    File expression_covariates
    File splicing_covariates
    File protein_covariates
    File? keep_samples

    # Gene and window selection
    Float qval_threshold = 0.05
    Int cis_window_bp = 1000000
    Int min_modalities = 2
    File? gene_list          # versioned ENSG IDs, one per line; overrides the qval filter
    Int? max_genes           # keep only the N most significant genes (useful for tests)
    Int min_shared_samples = 50

    # Optional gene picker. Set only_select_test_genes = true to run just SelectTestGenes
    # and write candidate genes for a first test; nothing else runs in that mode.
    Boolean only_select_test_genes = false
    Int test_min_junctions = 5
    Int test_max_junctions = 15
    Int test_min_variants = 8000
    Int test_max_variants = 14000
    Int test_n_candidates = 10

    # mvSuSiE settings (same defaults as the trans workflow)
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
    Int fit_memory_gb = 16

    # One image for every task (build it from cis-mvsusie-docker.zip)
    String docker_image
  }

  # 0. Optional: list candidate genes for a first test run, then stop.
  if (only_select_test_genes) {
    call SelectTestGenes {
      input:
        expression_cis_results = expression_cis_results,
        splicing_cis_results = splicing_cis_results,
        protein_cis_results = protein_cis_results,
        qval_threshold = qval_threshold,
        cis_window_bp = cis_window_bp,
        min_junctions = test_min_junctions,
        max_junctions = test_max_junctions,
        min_variants = test_min_variants,
        max_variants = test_max_variants,
        n_candidates = test_n_candidates,
        docker_image = docker_image
    }
  }

  if (!only_select_test_genes) {
    # 1. Sample IDs: everything is read as text and compared before any compute is spent.
    call CheckSamples {
      input:
        genome_dosage = genome_dosage,
        expression_phenotypes = expression_phenotypes,
        splicing_phenotypes = splicing_phenotypes,
        protein_phenotypes = protein_phenotypes,
        expression_covariates = expression_covariates,
        splicing_covariates = splicing_covariates,
        protein_covariates = protein_covariates,
        keep_samples = keep_samples,
        min_shared_samples = min_shared_samples,
        docker_image = docker_image
    }

    # 2. bgzip + tabix + lookup table for each phenotype BED (needed for fast per-gene reads)
    call PreparePhenotypeBed as IndexExpression {
      input: bed = expression_phenotypes, modality = "expression", docker_image = docker_image
    }
    call PreparePhenotypeBed as IndexSplicing {
      input: bed = splicing_phenotypes, modality = "splicing", docker_image = docker_image
    }
    call PreparePhenotypeBed as IndexProtein {
      input: bed = protein_phenotypes, modality = "protein", docker_image = docker_image
    }

    # 3. Choose genes and build the per-gene window, association and target tables
    call BuildCisWindows {
      input:
        expression_cis_results = expression_cis_results,
        splicing_cis_results = splicing_cis_results,
        protein_cis_results = protein_cis_results,
        expression_bed = IndexExpression.bed_bgz,
        splicing_bed = IndexSplicing.bed_bgz,
        protein_bed = IndexProtein.bed_bgz,
        shared_samples = CheckSamples.shared_samples,
        qval_threshold = qval_threshold,
        cis_window_bp = cis_window_bp,
        min_modalities = min_modalities,
        gene_list = gene_list,
        max_genes = max_genes,
        docker_image = docker_image
    }

    # 4. One pass per gene
    scatter (window_id in BuildCisWindows.window_ids) {
      call PrepareWindowGenotypes {
        input:
          window_id = window_id,
          trans_window_associations = BuildCisWindows.cis_window_associations,
          genome_dosage = genome_dosage,
          genome_dosage_tbi = genome_dosage_tbi,
          docker_image = docker_image
      }

      call PrepareWindowPhenotypes {
        input:
          window_id = window_id,
          trans_window_associations = BuildCisWindows.cis_window_associations,
          target_phenotypes = BuildCisWindows.target_phenotypes,
          expression_phenotypes = IndexExpression.bed_bgz,
          expression_phenotypes_tbi = IndexExpression.bed_tbi,
          expression_phenotype_lookup = IndexExpression.lookup,
          splicing_phenotypes = IndexSplicing.bed_bgz,
          splicing_phenotypes_tbi = IndexSplicing.bed_tbi,
          splicing_phenotype_lookup = IndexSplicing.lookup,
          protein_phenotypes = IndexProtein.bed_bgz,
          protein_phenotypes_tbi = IndexProtein.bed_tbi,
          protein_phenotype_lookup = IndexProtein.lookup,
          docker_image = docker_image
      }

      call PrepareMvSusieInput {
        input:
          window_id = window_id,
          window_manifest = PrepareWindowGenotypes.window_manifest,
          window_phenotypes_tsv = PrepareWindowPhenotypes.window_phenotypes,
          dosage = PrepareWindowGenotypes.window_dosage,
          phenotype_data = PrepareWindowPhenotypes.phenotype_data,
          expression_covariates = expression_covariates,
          splicing_covariates = splicing_covariates,
          protein_covariates = protein_covariates,
          keep_samples = CheckSamples.shared_samples,
          min_genotype_variance = min_genotype_variance,
          min_phenotype_variance = min_phenotype_variance,
          docker_image = docker_image
      }

      call FitMvSusie {
        input:
          window_id = window_id,
          prepared_window = PrepareMvSusieInput.prepared_window,
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
          memory_gb = fit_memory_gb,
          docker_image = docker_image
      }

      call SummarizeMvSusie {
        input:
          window_id = window_id,
          prepared_window = FitMvSusie.prepared_window_output,
          mvsusie_fit = FitMvSusie.mvsusie_fit,
          docker_image = docker_image
      }

      call PlotMvSusie {
        input:
          window_id = window_id,
          prepared_window = FitMvSusie.prepared_window_output,
          mvsusie_fit = FitMvSusie.mvsusie_fit,
          docker_image = docker_image
      }
    }

    # 5. Merge the per-gene tables (adds a window_id column when the table lacks one)
    call GatherTables as GatherVariantPip {
      input:
        files = SummarizeMvSusie.variant_pip,
        window_ids = BuildCisWindows.window_ids,
        output_name = "cis_mvsusie.variant_pip.tsv.gz",
        docker_image = docker_image
    }
    call GatherTables as GatherCredibleSets {
      input:
        files = SummarizeMvSusie.credible_sets,
        window_ids = BuildCisWindows.window_ids,
        output_name = "cis_mvsusie.credible_sets.tsv.gz",
        docker_image = docker_image
    }
    call GatherTables as GatherCredibleSetMembers {
      input:
        files = SummarizeMvSusie.credible_set_members,
        window_ids = BuildCisWindows.window_ids,
        output_name = "cis_mvsusie.credible_set_members.tsv.gz",
        docker_image = docker_image
    }
    call GatherTables as GatherComponentSupport {
      input:
        files = SummarizeMvSusie.component_feature_support,
        window_ids = BuildCisWindows.window_ids,
        output_name = "cis_mvsusie.component_feature_support.tsv.gz",
        docker_image = docker_image
    }
    call GatherTables as GatherWindowQc {
      input:
        files = SummarizeMvSusie.window_qc,
        window_ids = BuildCisWindows.window_ids,
        output_name = "cis_mvsusie.window_qc.tsv.gz",
        docker_image = docker_image
    }
  }

  output {
    File? test_candidates_all_modalities = SelectTestGenes.candidates_all_modalities
    File? test_candidates_no_protein = SelectTestGenes.candidates_no_protein
    File? suggested_gene_all_modalities = SelectTestGenes.suggested_gene_all_modalities
    File? suggested_gene_no_protein = SelectTestGenes.suggested_gene_no_protein

    File? sample_check = CheckSamples.sample_check
    File? cis_window_summary = BuildCisWindows.cis_window_summary
    Array[String]? window_ids = BuildCisWindows.window_ids

    File? merged_variant_pip = GatherVariantPip.merged
    File? merged_credible_sets = GatherCredibleSets.merged
    File? merged_credible_set_members = GatherCredibleSetMembers.merged
    File? merged_component_feature_support = GatherComponentSupport.merged
    File? merged_window_qc = GatherWindowQc.merged

    Array[File]? window_phenotypes = PrepareWindowPhenotypes.window_phenotypes
    Array[File]? window_prep_qc = PrepareWindowPhenotypes.window_qc
    Array[File]? prepared_windows = FitMvSusie.prepared_window_output
    Array[File]? mvsusie_fits = FitMvSusie.mvsusie_fit
    Array[File]? mashr_training = FitMvSusie.mashr_training
    Array[File]? greedy_L_history = FitMvSusie.greedy_L_history
    Array[File]? covariate_provenance = FitMvSusie.covariate_provenance
    Array[File]? run_stdout = FitMvSusie.run_stdout
    Array[File]? run_stderr = FitMvSusie.run_stderr
    Array[File]? session_info = FitMvSusie.session_info
    Array[File]? effect_plot_png = PlotMvSusie.effect_plot_png
    Array[File]? effect_plot_pdf = PlotMvSusie.effect_plot_pdf
    Array[File]? effect_plot_rds = PlotMvSusie.effect_plot_rds
  }
}

# ---------------------------------------------------------------------------
# SelectTestGenes (optional): suggests genes for a first test run, using only the cis
# permutation results. Two lists are written:
#   all_modalities: significant in expression, splicing and protein, with a moderate
#       junction count and variant count, not clipped at the chromosome start, ranked by
#       how close the three lead variants are to each other (then by closeness of
#       num_var to the median).
#   no_protein: significant in expression and splicing, with no protein result. Use it
#       to test the path where a modality is absent.
# The TSS is recovered as lead variant position minus start_distance.
# ---------------------------------------------------------------------------
task SelectTestGenes {
  input {
    File expression_cis_results
    File splicing_cis_results
    File protein_cis_results
    Float qval_threshold
    Int cis_window_bp
    Int min_junctions
    Int max_junctions
    Int min_variants
    Int max_variants
    Int n_candidates
    String docker_image
  }

  command <<<
    set -euo pipefail

    cat > select_test_genes.R <<'RSCRIPT'
    suppressPackageStartupMessages(library(data.table))
    log_msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))

    a <- commandArgs(trailingOnly = TRUE)
    modalities <- c("expression", "splicing", "protein")
    cis_paths <- setNames(a[1:3], modalities)
    qval_threshold <- as.numeric(a[4])
    window_bp <- as.integer(a[5])
    min_junctions <- as.integer(a[6])
    max_junctions <- as.integer(a[7])
    min_variants <- as.integer(a[8])
    max_variants <- as.integer(a[9])
    n_out <- as.integer(a[10])
    out_dir <- a[11]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    gene_from_feature <- function(mod, id) {
      if (mod == "expression") return(id)
      if (mod == "splicing") return(sub("^.*:", "", id))
      sub("^[^_]*_", "", id)
    }

    read_lead <- function(path, mod) {
      dt <- fread(
        cmd = sprintf("gzip -cdf %s", shQuote(path)),
        sep = "\t",
        colClasses = list(character = c("phenotype_id", "variant_id"))
      )
      need <- c("phenotype_id", "variant_id", "start_distance", "num_var", "qval")
      missing <- setdiff(need, names(dt))
      if (length(missing)) {
        stop(mod, " cis results are missing columns: ", paste(missing, collapse = ", "))
      }
      if (mod == "splicing" && "group_id" %in% names(dt)) {
        dt[, gene_id := as.character(group_id)]
      } else {
        dt[, gene_id := gene_from_feature(mod, phenotype_id)]
      }
      if (!"group_size" %in% names(dt)) dt[, group_size := NA_integer_]
      dt[, qval := suppressWarnings(as.numeric(qval))]
      dt[is.na(qval), qval := 1]
      setorder(dt, gene_id, qval)
      lead <- dt[!duplicated(gene_id)]
      lead[, lead_pos := as.integer(sub("^[^:]*:([0-9]+)_.*$", "\\1", variant_id))]
      lead[, chrom := sub(":.*$", "", variant_id)]
      lead[, tss := lead_pos - as.integer(start_distance)]
      out <- lead[, .(gene_id, qval, variant_id, lead_pos, chrom, tss,
                      num_var = as.integer(num_var), group_size = as.integer(group_size))]
      log_msg("%s: %d genes, %d with qval < %g.", mod, nrow(out), sum(out$qval < qval_threshold),
              qval_threshold)
      out
    }

    leads <- setNames(lapply(modalities, function(mod) read_lead(cis_paths[[mod]], mod)), modalities)
    e <- leads$expression[, .(gene_id, chrom, tss, num_var,
                              q_expression = qval, lead_expression = variant_id, pos_expression = lead_pos)]
    s <- leads$splicing[, .(gene_id, n_junctions = group_size,
                            q_splicing = qval, lead_splicing = variant_id, pos_splicing = lead_pos)]
    p <- leads$protein[, .(gene_id, q_protein = qval, lead_protein = variant_id, pos_protein = lead_pos)]
    w <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), list(e, s, p))

    junction_filter_available <- any(!is.na(w$n_junctions))
    if (!junction_filter_available) {
      log_msg("WARNING: no group_size values in the splicing results; the junction filter is skipped.")
    }
    common_filter <- function(d) {
      d <- d[!is.na(num_var) & num_var >= min_variants & num_var <= max_variants]
      d <- d[!is.na(tss) & tss > window_bp + 1L]
      if (junction_filter_available) {
        d <- d[!is.na(n_junctions) & n_junctions >= min_junctions & n_junctions <= max_junctions]
      }
      d
    }

    # Candidates with all three modalities significant
    A <- w[!is.na(q_expression) & !is.na(q_splicing) & !is.na(q_protein)]
    A <- A[q_expression < qval_threshold & q_splicing < qval_threshold & q_protein < qval_threshold]
    log_msg("%d genes are significant in all three modalities.", nrow(A))
    A <- common_filter(A)
    log_msg("%d of those pass the junction, variant and chromosome-start filters.", nrow(A))
    med_var <- median(w$num_var, na.rm = TRUE)
    if (nrow(A)) {
      A[, lead_span_bp := pmax(pos_expression, pos_splicing, pos_protein) -
          pmin(pos_expression, pos_splicing, pos_protein)]
      A[, n_distinct_leads := uniqueN(c(lead_expression, lead_splicing, lead_protein)), by = gene_id]
      A[, var_dist := abs(num_var - med_var)]
      setorder(A, lead_span_bp, var_dist, gene_id)
      A[, var_dist := NULL]
    }
    A <- head(A, n_out)

    # Candidates with no protein feature
    B <- w[!is.na(q_expression) & !is.na(q_splicing) & is.na(q_protein)]
    B <- B[q_expression < qval_threshold & q_splicing < qval_threshold]
    log_msg("%d genes are significant in expression and splicing with no protein result.", nrow(B))
    B <- common_filter(B)
    log_msg("%d of those pass the filters.", nrow(B))
    if (nrow(B)) {
      B[, worst_q := pmax(q_expression, q_splicing)]
      B[, lead_span_bp := abs(pos_expression - pos_splicing)]
      setorder(B, lead_span_bp, worst_q, gene_id)
      B[, worst_q := NULL]
    }
    B <- head(B, n_out)

    cols_all <- c("gene_id", "chrom", "tss", "num_var", "n_junctions",
                  "q_expression", "q_splicing", "q_protein",
                  "lead_expression", "lead_splicing", "lead_protein",
                  "lead_span_bp", "n_distinct_leads")
    cols_np <- c("gene_id", "chrom", "tss", "num_var", "n_junctions",
                 "q_expression", "q_splicing",
                 "lead_expression", "lead_splicing", "lead_span_bp")
    write_table <- function(d, cols, path) {
      if (!nrow(d)) {
        d <- as.data.table(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols))
      }
      fwrite(d[, ..cols], path, sep = "\t", quote = FALSE)
    }
    write_table(A, cols_all, file.path(out_dir, "test_candidates_all_modalities.tsv"))
    write_table(B, cols_np, file.path(out_dir, "test_candidates_no_protein.tsv"))
    writeLines(if (nrow(A)) A$gene_id[1L] else character(), file.path(out_dir, "suggested_gene_all_modalities.txt"))
    writeLines(if (nrow(B)) B$gene_id[1L] else character(), file.path(out_dir, "suggested_gene_no_protein.txt"))
    if (!nrow(A)) log_msg("WARNING: no all-modality candidates. Loosen the junction or variant limits.")
    if (!nrow(B)) log_msg("WARNING: no no-protein candidates. Loosen the junction or variant limits.")
    RSCRIPT

    mkdir -p out
    Rscript select_test_genes.R \
      '~{expression_cis_results}' '~{splicing_cis_results}' '~{protein_cis_results}' \
      '~{qval_threshold}' '~{cis_window_bp}' \
      '~{min_junctions}' '~{max_junctions}' '~{min_variants}' '~{max_variants}' \
      '~{n_candidates}' out
    echo "All-modality suggestion: $(cat out/suggested_gene_all_modalities.txt)" >&2
    echo "No-protein suggestion: $(cat out/suggested_gene_no_protein.txt)" >&2
    echo "Pass a suggested_gene_*.txt output file as gene_list and set only_select_test_genes = false." >&2
  >>>

  output {
    File candidates_all_modalities = "out/test_candidates_all_modalities.tsv"
    File candidates_no_protein = "out/test_candidates_no_protein.tsv"
    File suggested_gene_all_modalities = "out/suggested_gene_all_modalities.txt"
    File suggested_gene_no_protein = "out/suggested_gene_no_protein.txt"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "8 GiB"
    disks: "local-disk 50 SSD"
    bootDiskSizeGb: 30
  }
}

# ---------------------------------------------------------------------------
# CheckSamples: compare sample IDs across every input. IDs are always treated as text,
# whitespace is trimmed, and a leading X before a digit is removed (the same
# normalization the R scripts use). Fails if too few samples are shared.
# ---------------------------------------------------------------------------
task CheckSamples {
  input {
    File genome_dosage
    File expression_phenotypes
    File splicing_phenotypes
    File protein_phenotypes
    File expression_covariates
    File splicing_covariates
    File protein_covariates
    File? keep_samples
    Int min_shared_samples
    String docker_image
  }

  command <<<
    set -euo pipefail
    export LC_ALL=C
    log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2; }
    mkdir -p ids

    first_line() {
      { gzip -cdf "$1" || true; } | head -n 1 | tr -d '\r'
    }

    # extract_ids <file> <first sample column> <name>
    extract_ids() {
      first_line "$1" | awk -F '\t' -v start="$2" '
        {
          for (i = start; i <= NF; i++) {
            id = $i
            gsub(/^[ \t]+|[ \t]+$/, "", id)
            if (id ~ /^X[0-9]/) id = substr(id, 2)
            print id
          }
        }' > "ids/$3.raw"
      if [[ -n "$(sort "ids/$3.raw" | uniq -d)" ]]; then
        log "ERROR: duplicate sample IDs after normalization in $3."
        exit 1
      fi
      sort -u "ids/$3.raw" > "ids/$3.txt"
      rm "ids/$3.raw"
    }

    log "Reading sample IDs from every header."
    extract_ids '~{genome_dosage}' 5 dosage
    extract_ids '~{expression_phenotypes}' 5 expression_bed
    extract_ids '~{splicing_phenotypes}' 5 splicing_bed
    extract_ids '~{protein_phenotypes}' 5 protein_bed
    extract_ids '~{expression_covariates}' 2 expression_cov
    extract_ids '~{splicing_covariates}' 2 splicing_cov
    extract_ids '~{protein_covariates}' 2 protein_cov

    keep='~{default="" keep_samples}'
    if [[ -n "$keep" ]]; then
      tr -d '\r' < "$keep" | awk '
        {
          gsub(/^[ \t]+|[ \t]+$/, "")
          if ($0 ~ /^X[0-9]/) $0 = substr($0, 2)
          if ($0 != "" && tolower($0) != "research_id") print
        }' | sort -u > ids/keep.txt
    fi

    cp ids/dosage.txt shared.txt
    for f in ids/*.txt; do
      comm -12 shared.txt "$f" > shared.next
      mv shared.next shared.txt
    done
    n_shared="$(wc -l < shared.txt)"
    printf '%s\n' "$n_shared" > n_shared.txt

    {
      printf 'source\tn_samples\tn_in_shared\n'
      for f in ids/*.txt; do
        name="$(basename "$f" .txt)"
        printf '%s\t%s\t%s\n' "$name" "$(wc -l < "$f")" "$(comm -12 "$f" shared.txt | wc -l)"
      done
      printf 'shared\t%s\t%s\n' "$n_shared" "$n_shared"
    } > sample_check.tsv
    cat sample_check.tsv >&2

    if (( n_shared < ~{min_shared_samples} )); then
      log "ERROR: only ${n_shared} shared samples; at least ~{min_shared_samples} are required."
      exit 1
    fi
    cp shared.txt shared_samples.txt
    log "Sample check complete: ${n_shared} shared samples."
  >>>

  output {
    File sample_check = "sample_check.tsv"
    File shared_samples = "shared_samples.txt"
    Int n_shared_samples = read_int("n_shared.txt")
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: 30
  }
}

# ---------------------------------------------------------------------------
# PreparePhenotypeBed: re-bgzip (sorted), tabix-index and build the lookup table
# (phenotype_id, chrom, start, end) that prepare_trans_window.R uses for indexed access.
# ---------------------------------------------------------------------------
task PreparePhenotypeBed {
  input {
    File bed
    String modality
    String docker_image
  }

  command <<<
    set -euo pipefail
    export LC_ALL=C
    log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2; }

    log "Sorting, bgzipping and indexing the ~{modality} BED."
    { gzip -cdf '~{bed}' || true; } | tr -d '\r' | {
      IFS= read -r header
      if [[ "${header}" != "#"* ]]; then
        echo "ERROR: the first line of the ~{modality} BED must start with '#'." >&2
        exit 1
      fi
      printf '%s\n' "${header}"
      sort -t $'\t' -k1,1 -k2,2n -k3,3n -S 4G -T .
    } | bgzip -c > '~{modality}.bed.gz'
    tabix -f -p bed '~{modality}.bed.gz'
    test -s '~{modality}.bed.gz'
    test -s '~{modality}.bed.gz.tbi'

    log "Building the lookup table."
    printf 'phenotype_id\tchrom\tstart\tend\n' > '~{modality}.lookup.tsv'
    gzip -cd '~{modality}.bed.gz' \
      | awk -F '\t' 'BEGIN { OFS = "\t" } NR > 1 { print $4, $1, $2, $3 }' \
      >> '~{modality}.lookup.tsv'
    n_rows="$(( $(wc -l < '~{modality}.lookup.tsv') - 1 ))"
    n_unique="$(tail -n +2 '~{modality}.lookup.tsv' | cut -f1 | sort -u | wc -l)"
    if (( n_rows < 1 )) || (( n_rows != n_unique )); then
      log "ERROR: ~{modality} BED has ${n_rows} rows but ${n_unique} unique phenotype IDs."
      exit 1
    fi
    log "Indexed ${n_rows} ~{modality} features."
  >>>

  output {
    File bed_bgz = modality + ".bed.gz"
    File bed_tbi = modality + ".bed.gz.tbi"
    File lookup = modality + ".lookup.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 200 SSD"
    bootDiskSizeGb: 30
  }
}

# ---------------------------------------------------------------------------
# BuildCisWindows: pick genes and write the tables the reused trans scripts expect.
#   cis_window_associations.tsv.gz: window_id chrom start end modality molecular_trait_id p_value
#   target_phenotypes.tsv:          window_id modality phenotype_id (expression + splicing only)
# Gene IDs are taken from: expression -> the ID; splicing -> cis group_id / last ":" field;
# protein -> the ID after the UniProt prefix. The TSS is the BED "end" column (1-based).
# Window = [TSS - W, TSS + W] in 1-based positions, stored 0-based half-open as the
# manifest convention requires: start = max(0, TSS - W - 1), end = TSS + W.
# ---------------------------------------------------------------------------
task BuildCisWindows {
  input {
    File expression_cis_results
    File splicing_cis_results
    File protein_cis_results
    File expression_bed
    File splicing_bed
    File protein_bed
    File shared_samples
    Float qval_threshold
    Int cis_window_bp
    Int min_modalities
    File? gene_list
    Int? max_genes
    String docker_image
  }

  command <<<
    set -euo pipefail
    echo "[$(date -u +'%Y-%m-%dT%H:%M:%SZ')] Building cis windows. Shared samples: $(wc -l < '~{shared_samples}')." >&2

    cat > build_cis_windows.R <<'RSCRIPT'
    suppressPackageStartupMessages(library(data.table))
    log_msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))

    a <- commandArgs(trailingOnly = TRUE)
    modalities <- c("expression", "splicing", "protein")
    cis_paths <- setNames(a[1:3], modalities)
    bed_paths <- setNames(a[4:6], modalities)
    qval_threshold <- as.numeric(a[7])
    window_bp <- as.integer(a[8])
    min_modalities <- as.integer(a[9])
    gene_list_path <- a[10]
    max_genes <- as.integer(a[11])
    out_dir <- a[12]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    gene_from_feature <- function(mod, id) {
      if (mod == "expression") return(id)
      if (mod == "splicing") return(sub("^.*:", "", id))
      sub("^[^_]*_", "", id)
    }

    read_cis <- function(path, mod) {
      dt <- fread(
        cmd = sprintf("gzip -cdf %s", shQuote(path)),
        sep = "\t",
        colClasses = list(character = "phenotype_id")
      )
      need <- c("phenotype_id", "pval_nominal", "qval")
      missing <- setdiff(need, names(dt))
      if (length(missing)) {
        stop(mod, " cis results are missing columns: ", paste(missing, collapse = ", "))
      }
      dt[, phenotype_id := as.character(phenotype_id)]
      dt[, qval := suppressWarnings(as.numeric(qval))]
      dt[, pval_nominal := suppressWarnings(as.numeric(pval_nominal))]
      dt[is.na(qval), qval := 1]
      dt[is.na(pval_nominal), pval_nominal := 1]
      if (mod == "splicing" && "group_id" %in% names(dt)) {
        dt[, gene_id := as.character(group_id)]
      } else {
        dt[, gene_id := gene_from_feature(mod, phenotype_id)]
      }
      log_msg("%s cis results: %d rows, %d genes, %d with qval < %g.",
              mod, nrow(dt), uniqueN(dt$gene_id), uniqueN(dt[qval < qval_threshold, gene_id]),
              qval_threshold)
      dt[, .(phenotype_id, gene_id, pval_nominal, qval)]
    }

    read_bed_features <- function(path, mod) {
      dt <- fread(
        cmd = sprintf("gzip -cdf %s | cut -f1-4", shQuote(path)),
        sep = "\t", header = TRUE,
        colClasses = c("character", "integer", "integer", "character")
      )
      setnames(dt, c("chrom", "start", "end", "feature_id"))
      if (anyDuplicated(dt$feature_id)) stop(mod, " BED has duplicate feature IDs.")
      dt[, gene_id := gene_from_feature(mod, feature_id)]
      dt[, modality := mod]
      log_msg("%s BED: %d features, %d genes.", mod, nrow(dt), uniqueN(dt$gene_id))
      dt[]
    }

    cis <- setNames(lapply(modalities, function(mod) read_cis(cis_paths[[mod]], mod)), modalities)
    gene_cis <- rbindlist(lapply(modalities, function(mod) {
      cis[[mod]][, .(min_qval = min(qval), best_p = min(pval_nominal)), by = gene_id][, modality := mod]
    }))
    overall <- gene_cis[, .(min_qval = min(min_qval)), by = gene_id]

    if (gene_list_path != "NONE") {
      genes <- unique(trimws(readLines(gene_list_path, warn = FALSE)))
      genes <- genes[nzchar(genes)]
      log_msg("Using %d genes from the supplied list; the qval filter is ignored.", length(genes))
    } else {
      genes <- overall[min_qval < qval_threshold, gene_id]
      log_msg("%d genes have qval < %g in at least one modality.", length(genes), qval_threshold)
    }

    features <- rbindlist(lapply(modalities, function(mod) read_bed_features(bed_paths[[mod]], mod)))
    n_by <- features[gene_id %in% genes, .N, by = .(gene_id, modality)]
    gene_tbl <- data.table(gene_id = genes)
    for (mod in modalities) {
      n_mod <- n_by[modality == mod]
      idx <- match(gene_tbl$gene_id, n_mod$gene_id)
      set(gene_tbl, j = paste0("n_", mod), value = ifelse(is.na(idx), 0L, n_mod$N[idx]))
    }
    gene_tbl[, n_present := (n_expression > 0L) + (n_splicing > 0L) + (n_protein > 0L)]
    keep_flag <- gene_tbl$n_present >= min_modalities &
      (gene_tbl$n_expression > 0L | gene_tbl$n_splicing > 0L)
    log_msg("%d of %d genes have features in >= %d modalities (including expression or splicing).",
            sum(keep_flag), nrow(gene_tbl), min_modalities)
    gene_tbl <- gene_tbl[keep_flag]
    gene_tbl <- merge(gene_tbl, overall, by = "gene_id", all.x = TRUE, sort = FALSE)
    setorder(gene_tbl, min_qval, gene_id, na.last = TRUE)
    if (!is.na(max_genes) && max_genes > 0L) {
      gene_tbl <- head(gene_tbl, max_genes)
      log_msg("Keeping the %d most significant genes.", nrow(gene_tbl))
    }
    if (!nrow(gene_tbl)) stop("No genes passed the selection.")
    if (!all(grepl("^[A-Za-z0-9._-]+$", gene_tbl$gene_id))) {
      stop("Gene IDs may contain letters, numbers, periods, underscores and hyphens only.")
    }

    feat <- features[gene_id %in% gene_tbl$gene_id]
    feat[, prio := match(modality, modalities)]
    setorder(feat, gene_id, prio)
    pos <- feat[, .(
      window_chrom = chrom[1L],
      tss = end[1L],
      n_distinct_tss = uniqueN(paste(chrom, end))
    ), by = gene_id]
    if (any(pos$n_distinct_tss > 1L)) {
      log_msg("WARNING: %d genes have different TSS positions across features; using the first modality that has the gene (expression, splicing, protein).",
              sum(pos$n_distinct_tss > 1L))
    }
    pos[, window_start := pmax(0L, tss - window_bp - 1L)]
    pos[, window_end := tss + window_bp]

    cis_feat <- unique(rbindlist(lapply(modalities, function(mod) {
      cis[[mod]][, .(modality = mod, feature_id = phenotype_id, own_p = pval_nominal)]
    })), by = c("modality", "feature_id"))
    feat <- merge(feat, pos[, .(gene_id, window_chrom, window_start, window_end)],
                  by = "gene_id", sort = FALSE)
    feat <- merge(feat, cis_feat, by = c("modality", "feature_id"), all.x = TRUE, sort = FALSE)
    feat <- merge(feat, gene_cis[, .(gene_id, modality, best_p)],
                  by = c("gene_id", "modality"), all.x = TRUE, sort = FALSE)
    feat[, p_value := pmin(pmax(fcoalesce(own_p, best_p, 1), 0), 1)]
    feat[, gene_rank := match(gene_id, gene_tbl$gene_id)]
    setorder(feat, gene_rank, prio, feature_id)

    assoc <- feat[, .(
      window_id = gene_id,
      chrom = window_chrom,
      start = as.integer(window_start),
      end = as.integer(window_end),
      modality = modality,
      molecular_trait_id = feature_id,
      p_value = p_value
    )]
    targets <- assoc[modality %in% c("expression", "splicing"),
                     .(window_id, modality, phenotype_id = molecular_trait_id)]

    qwide <- dcast(gene_cis, gene_id ~ modality, value.var = "min_qval")
    for (mod in modalities) {
      if (!mod %in% names(qwide)) set(qwide, j = mod, value = NA_real_)
    }
    setnames(qwide, modalities, paste0("qval_", modalities))
    summary_dt <- merge(gene_tbl[, .(gene_id, n_expression, n_splicing, n_protein)],
                        pos[, .(gene_id, chrom = window_chrom, tss, window_start, window_end, n_distinct_tss)],
                        by = "gene_id", sort = FALSE)
    summary_dt <- merge(summary_dt, qwide, by = "gene_id", all.x = TRUE, sort = FALSE)
    summary_dt[, gene_rank := match(gene_id, gene_tbl$gene_id)]
    setorder(summary_dt, gene_rank)
    summary_dt[, gene_rank := NULL]
    setnames(summary_dt, "gene_id", "window_id")

    fwrite(assoc, file.path(out_dir, "cis_window_associations.tsv"), sep = "\t", quote = FALSE)
    fwrite(targets, file.path(out_dir, "target_phenotypes.tsv"), sep = "\t", quote = FALSE)
    fwrite(summary_dt, file.path(out_dir, "cis_window_summary.tsv"), sep = "\t", quote = FALSE)
    writeLines(gene_tbl$gene_id, file.path(out_dir, "window_ids.txt"))
    log_msg("Wrote %d windows, %d association rows and %d target rows.",
            nrow(gene_tbl), nrow(assoc), nrow(targets))
    RSCRIPT

    mkdir -p out
    Rscript build_cis_windows.R \
      '~{expression_cis_results}' '~{splicing_cis_results}' '~{protein_cis_results}' \
      '~{expression_bed}' '~{splicing_bed}' '~{protein_bed}' \
      '~{qval_threshold}' '~{cis_window_bp}' '~{min_modalities}' \
      '~{default="NONE" gene_list}' '~{default=0 max_genes}' out
    gzip -f out/cis_window_associations.tsv
    test -s out/cis_window_associations.tsv.gz
    test -s out/target_phenotypes.tsv
    test -s out/window_ids.txt
    test -s out/cis_window_summary.tsv
  >>>

  output {
    File cis_window_associations = "out/cis_window_associations.tsv.gz"
    File target_phenotypes = "out/target_phenotypes.tsv"
    File cis_window_summary = "out/cis_window_summary.tsv"
    Array[String] window_ids = read_lines("out/window_ids.txt")
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: 30
  }
}

# ---------------------------------------------------------------------------
# PrepareWindowGenotypes: unchanged from the collaborator's workflow, except that the
# image is a parameter.
# ---------------------------------------------------------------------------
task PrepareWindowGenotypes {
  input {
    File trans_window_associations
    String window_id
    File genome_dosage
    File genome_dosage_tbi
    String docker_image
  }

  command <<<
    set -euo pipefail

    output_prefix='~{window_id}'
    if [[ ! "$output_prefix" =~ ^[A-Za-z0-9._-]+$ ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: window_id contains unsafe filename characters." >&2
      exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting genotype preparation for ~{window_id}."
    mkdir -p output

    dosage_name="window_input.dose.tsv.gz"
    ln -sf "~{genome_dosage}" "${dosage_name}"
    ln -sf "~{genome_dosage_tbi}" "${dosage_name}.tbi"
    test -s "${dosage_name}"
    test -s "${dosage_name}.tbi"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Local dosage: ${dosage_name}."
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Local index: ${dosage_name}.tbi."

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Resolving the locus coordinates."
    window_row="$(awk -F '\t' -v requested_id='~{window_id}' '
      NR == 1 {
        for (i = 1; i <= NF; i++) column[$i] = i
        next
      }
      $(column["window_id"]) == requested_id {
        row = $(column["chrom"]) "\t" $(column["start"]) "\t" $(column["end"])
        if (!seen[row]++) {
          selected = row
          matches++
        }
      }
      END {
        if (matches != 1) exit 1
        print selected
      }
    ' <(gzip -cd ~{trans_window_associations}))"

    IFS=$'\t' read -r window_chrom window_start window_end <<< "${window_row}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Interval: ${window_chrom}:$((window_start + 1))-${window_end}."
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Reading the dosage header."
    tabix -H "${dosage_name}" > dosage_header.raw.tsv
    if [[ -s dosage_header.raw.tsv ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using the header returned by tabix."
    else
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] tabix returned no header. Reading the first compressed line."
      set +o pipefail
      gzip -cd -- "${dosage_name}" \
        | awk 'NR == 1 {print; exit}' \
        > dosage_header.raw.tsv
      set -o pipefail
    fi
    test -s dosage_header.raw.tsv
    header_line="$(tail -n 1 dosage_header.raw.tsv)"
    header_line="${header_line#\#}"
    printf '%s\n' "${header_line}" > "output/${output_prefix}.window_dosage.tsv"
    if ! awk -F '\t' '
      NR == 1 {
        valid = ($1 == "CHROM" && $2 == "POS" && $3 == "REF" && $4 == "ALT")
        exit !valid
      }
    ' "output/${output_prefix}.window_dosage.tsv"; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Dosage header must start with CHROM, POS, REF, and ALT." >&2
      exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Dosage header validation complete."

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting all locus genotype rows."
    tabix "${dosage_name}" \
      "${window_chrom}:$((window_start + 1))-${window_end}" \
      >> "output/${output_prefix}.window_dosage.tsv"
    test "$(wc -l < "output/${output_prefix}.window_dosage.tsv")" -gt 1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracted variants: $(( $(wc -l < "output/${output_prefix}.window_dosage.tsv") - 1 ))."

    {
      printf 'window_id\tchrom\tstart\tend\tdosage_file\n'
      printf '%s\t%s\t%s\t%s\t%s\n' \
        '~{window_id}' \
        "${window_chrom}" \
        "${window_start}" \
        "${window_end}" \
        "${output_prefix}.window_dosage.tsv"
    } > "output/${output_prefix}.window_manifest.tsv"
    test -s "output/${output_prefix}.window_manifest.tsv"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Genotype preparation complete."
  >>>

  output {
    File window_dosage = "output/" + window_id + ".window_dosage.tsv"
    File window_manifest = "output/" + window_id + ".window_manifest.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
    bootDiskSizeGb: 30
  }
}

# ---------------------------------------------------------------------------
# PrepareWindowPhenotypes: the collaborator's script, called with indexed inputs for
# every modality and top-N caps high enough to keep every feature of the gene.
# ---------------------------------------------------------------------------
task PrepareWindowPhenotypes {
  input {
    String window_id
    File trans_window_associations
    File target_phenotypes
    File expression_phenotypes
    File expression_phenotypes_tbi
    File expression_phenotype_lookup
    File splicing_phenotypes
    File splicing_phenotypes_tbi
    File splicing_phenotype_lookup
    File protein_phenotypes
    File protein_phenotypes_tbi
    File protein_phenotype_lookup
    Int top_n = 100000
    String docker_image
  }

  command <<<
    set -euo pipefail

    output_prefix='~{window_id}'
    if [[ ! "$output_prefix" =~ ^[A-Za-z0-9._-]+$ ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: window_id contains unsafe filename characters." >&2
      exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting joint phenotype preparation for ~{window_id}."

    Rscript /opt/mvsusie/scripts/prepare_trans_window.R \
      --window-id ~{window_id} \
      --trans-associations ~{trans_window_associations} \
      --expression-phenotypes ~{expression_phenotypes} \
      --splicing-phenotypes ~{splicing_phenotypes} \
      --protein-phenotypes ~{protein_phenotypes} \
      --expression-phenotypes-tbi ~{expression_phenotypes_tbi} \
      --expression-phenotype-lookup ~{expression_phenotype_lookup} \
      --splicing-phenotypes-tbi ~{splicing_phenotypes_tbi} \
      --splicing-phenotype-lookup ~{splicing_phenotype_lookup} \
      --protein-phenotypes-tbi ~{protein_phenotypes_tbi} \
      --protein-phenotype-lookup ~{protein_phenotype_lookup} \
      --target-phenotypes ~{target_phenotypes} \
      --top-n-expression ~{top_n} \
      --top-n-splicing ~{top_n} \
      --top-n-protein ~{top_n} \
      --output-dir output
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Validating joint phenotype outputs."
    test -s "output/${output_prefix}.window_phenotypes.tsv"
    test -s "output/${output_prefix}.window_phenotypes.bed.gz"
    test -s "output/${output_prefix}.window_qc.tsv"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Joint phenotype preparation complete."
  >>>

  output {
    File window_phenotypes = "output/" + window_id + ".window_phenotypes.tsv"
    File phenotype_data = "output/" + window_id + ".window_phenotypes.bed.gz"
    File window_qc = "output/" + window_id + ".window_qc.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
    bootDiskSizeGb: 30
  }
}

# ---------------------------------------------------------------------------
# PrepareMvSusieInput, FitMvSusie, SummarizeMvSusie, PlotMvSusie: the collaborator's
# tasks. The only changes are the memory parameter on FitMvSusie and the log message
# that reports it.
# ---------------------------------------------------------------------------
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
    Int memory_gb
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
    log "Starting joint mvSuSiE for ~{window_id}. The task has ~{memory_gb} GiB of memory."
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
    memory: "~{memory_gb} GiB"
    disks: "local-disk 500 SSD"
    bootDiskSizeGb: 30
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
    bootDiskSizeGb: 30
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
    bootDiskSizeGb: 30
  }
}

# ---------------------------------------------------------------------------
# GatherTables: concatenate per-gene tables. Empty files are skipped. A window_id
# column is prepended unless the tables already have one.
# ---------------------------------------------------------------------------
task GatherTables {
  input {
    Array[File] files
    Array[String] window_ids
    String output_name
    String docker_image
  }

  command <<<
    set -euo pipefail
    files=(~{sep=" " files})
    ids=(~{sep=" " window_ids})
    if (( ${#files[@]} != ${#ids[@]} )); then
      echo "ERROR: ${#files[@]} files but ${#ids[@]} window IDs." >&2
      exit 1
    fi

    : > merged.tsv
    add_id=-1
    for i in "${!files[@]}"; do
      { gzip -cdf "${files[$i]}" || true; } > current.tsv
      [[ -s current.tsv ]] || continue
      if (( add_id < 0 )); then
        if awk -F '\t' 'NR == 1 { for (j = 1; j <= NF; j++) if ($j == "window_id") found = 1; exit } END { exit !found }' current.tsv; then
          add_id=0
        else
          add_id=1
          printf 'window_id\t' >> merged.tsv
        fi
        head -n 1 current.tsv >> merged.tsv
      fi
      if (( add_id == 1 )); then
        tail -n +2 current.tsv | awk -v id="${ids[$i]}" 'BEGIN { OFS = "\t" } { print id, $0 }' >> merged.tsv
      else
        tail -n +2 current.tsv >> merged.tsv
      fi
    done
    gzip -c merged.tsv > '~{output_name}'
    test -e '~{output_name}'
  >>>

  output {
    File merged = output_name
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 50 SSD"
    bootDiskSizeGb: 30
  }
}
