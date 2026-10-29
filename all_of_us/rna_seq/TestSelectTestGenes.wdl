version 1.0

workflow TestSelectTestGenes {
  input {
    File expression_cis_results
    File splicing_cis_results
    File protein_cis_results
    Float qval_threshold = 0.05
    Int cis_window_bp = 1000000
    Int min_junctions = 5
    Int max_junctions = 15
    Int min_variants = 8000
    Int max_variants = 14000
    Int n_candidates = 10
    String docker_image
  }

  call SelectTestGenes {
    input:
      expression_cis_results = expression_cis_results,
      splicing_cis_results = splicing_cis_results,
      protein_cis_results = protein_cis_results,
      qval_threshold = qval_threshold,
      cis_window_bp = cis_window_bp,
      min_junctions = min_junctions,
      max_junctions = max_junctions,
      min_variants = min_variants,
      max_variants = max_variants,
      n_candidates = n_candidates,
      docker_image = docker_image
  }

  output {
    File candidates_all_modalities = SelectTestGenes.candidates_all_modalities
    File candidates_no_protein = SelectTestGenes.candidates_no_protein
    File suggested_gene_all_modalities = SelectTestGenes.suggested_gene_all_modalities
    File suggested_gene_no_protein = SelectTestGenes.suggested_gene_no_protein
  }
}

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