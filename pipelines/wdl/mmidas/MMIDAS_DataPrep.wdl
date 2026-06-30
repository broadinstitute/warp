version 1.0

workflow MMIDAS_DataPrep {

  meta {
    description: "Stage 1 of the MMIDAS pipeline. Loads raw Mouse Smart-seq ALM/VISp count matrices from the Allen Brain Atlas, filters to neuronal cells, normalizes to log-CPM, and writes a single AnnData .h5ad file for use by MMIDAS_Train."
    allowNestedInputs: true
  }

  String pipeline_version = "1.0.0"

  input {
    # ── Raw Allen Brain Atlas Smart-seq files ────────────────────────────────
    File visp_exon_matrix   # mouse_VISp_2018-06-14_exon-matrix.csv
    File visp_samples       # mouse_VISp_2018-06-14_samples-columns.csv
    File alm_exon_matrix    # mouse_ALM_2018-06-14_exon-matrix.csv
    File alm_samples        # mouse_ALM_2018-06-14_samples-columns.csv
    File genes_rows         # mouse_ALM_2018-06-14_genes-rows.csv  (full gene list)
    File selected_genes     # genes_SS_ALM-VISp.csv  (selected gene subset ~1252 genes)

    # ── Output filename ──────────────────────────────────────────────────────
    String output_basename = "Mouse_ALM-VISp_cpm"

    # ── Optional filters ─────────────────────────────────────────────────────
    String remove_clusters  = "Low Quality,CR Lhx5,Meis2 Adamts19"
    String neuronal_classes = "GABAergic,Glutamatergic"

    # ── Runtime ──────────────────────────────────────────────────────────────
    String docker    = "us.gcr.io/broad-gotc-prod/mmidas:1.0.0-0.1.0-1782844522"
    Int    disk_size = 100
    Int    mem_size  = 32
    Int    cpu       = 4
  }

  call DataPrep {
    input:
      visp_exon_matrix  = visp_exon_matrix,
      visp_samples      = visp_samples,
      alm_exon_matrix   = alm_exon_matrix,
      alm_samples       = alm_samples,
      genes_rows        = genes_rows,
      selected_genes    = selected_genes,
      output_basename   = output_basename,
      remove_clusters   = remove_clusters,
      neuronal_classes  = neuronal_classes,
      docker            = docker,
      disk_size         = disk_size,
      mem_size          = mem_size,
      cpu               = cpu
  }

  output {
    File   preprocessed_h5ad     = DataPrep.preprocessed_h5ad
    String pipeline_version_out  = pipeline_version
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task: DataPrep
#
# Runs 01_data_prep.py which:
#   - Loads VISp and ALM count matrices
#   - Filters to neuronal classes (GABAergic + Glutamatergic by default)
#   - Removes low-quality / rare clusters
#   - Normalizes: log1p(counts / sum * 1e6)  (log-CPM)
#   - Subsets to the selected gene list
#   - Writes output as AnnData .h5ad
# ──────────────────────────────────────────────────────────────────────────────
task DataPrep {
  input {
    File   visp_exon_matrix
    File   visp_samples
    File   alm_exon_matrix
    File   alm_samples
    File   genes_rows
    File   selected_genes
    String output_basename
    String remove_clusters
    String neuronal_classes
    String docker
    Int    disk_size
    Int    mem_size
    Int    cpu
  }

  parameter_meta {
    visp_exon_matrix: "Raw VISp exon count matrix CSV (genes × cells)."
    visp_samples:     "VISp cell metadata CSV (samples-columns)."
    alm_exon_matrix:  "Raw ALM exon count matrix CSV (genes × cells)."
    alm_samples:      "ALM cell metadata CSV (samples-columns)."
    genes_rows:       "Full gene list CSV (genes-rows) matching the count matrices."
    selected_genes:   "Selected gene subset CSV (~1252 genes for Smart-seq ALM/VISp)."
    output_basename:  "Basename for the output .h5ad file (no extension)."
    remove_clusters:  "Comma-separated cluster names to exclude (e.g. 'Low Quality,CR Lhx5')."
    neuronal_classes: "Comma-separated cell classes to retain (default: GABAergic,Glutamatergic)."
    docker:           "MMIDAS Docker image URI."
    disk_size:        "Boot disk size in GB."
    mem_size:         "Memory in GB."
    cpu:              "Number of CPU cores."
  }

  command <<<
    set -euo pipefail

    python3 /usr/local/01_data_prep.py \
      --visp_exon_matrix  "~{visp_exon_matrix}" \
      --visp_samples      "~{visp_samples}" \
      --alm_exon_matrix   "~{alm_exon_matrix}" \
      --alm_samples       "~{alm_samples}" \
      --genes_rows        "~{genes_rows}" \
      --selected_genes    "~{selected_genes}" \
      --output_h5ad       "~{output_basename}.h5ad" \
      --remove_clusters   "~{remove_clusters}" \
      --neuronal_classes  "~{neuronal_classes}"
  >>>

  output {
    File preprocessed_h5ad = "~{output_basename}.h5ad"
  }

  runtime {
    docker:  docker
    disks:   "local-disk ~{disk_size} HDD"
    memory:  "~{mem_size} GiB"
    cpu:     cpu
    preemptible: 3
  }
}
