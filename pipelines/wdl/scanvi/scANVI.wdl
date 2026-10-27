version 1.0

import "../../../tasks/wdl/Utilities.wdl" as utils

workflow scANVI {

  meta {
    description: "Pipeline for cell type label transfer using SCVI and SCANVI models. Labels one or more query files against a single annotated reference: each query's single-cell RNA (GEX) and (optionally) ATAC data are integrated with the reference to transfer cell type labels via semi-supervised deep generative models, independently of the other queries and with their own outputs. When no ATAC h5ads are provided, the pipeline runs in GEX-only mode and trains/annotates from the reference atlas using GEX and reference alone."
    allowNestedInputs: true
  }

  input {
      # Query files, paired by position: input_ids[i] prefixes every output of gex_h5ads[i]
      # (and atac_h5ads[i]). A single query is a one-element array.
      Array[String] input_ids
      Array[File] gex_h5ads

      # Annotated reference h5ad, shared by every query
      File ref_h5ad

      # ATAC cell-by-bin h5ads for multiome mode: one per gex_h5ads entry, in the same order.
      # Empty (default) runs every query in GEX-only mode.
      Array[File] atac_h5ads = []

      # Optional cap on SCVI/SCANVI training epochs (applies to both multiome and
      # GEX-only modes). When unset, the container default (500) is used.
      Int? max_epochs

      # SCVI/SCANVI minibatch size (SGD). Default 128 (scvi-tools' own default, reproducing prior
      # behavior). Lower it (e.g. 8) to fit a high-cardinality reference on a small GPU; raise it on
      # a large-VRAM cloud GPU. Activation memory scales with (batch_size x number of labels).
      Int batch_size = 128

      # Optional random seed for SCVI/SCANVI training. Unset (default) keeps training
      # nondeterministic; set it to make a run reproducible (tests pin one).
      Int? seed

      # Reference adaptation. When the reference is an AIT-schema atlas, these select
      # which obs columns become the cell-type label and the batch. When unset, AIT
      # references default to subclass/donor_id and PBMC-style references to
      # final_annotation/batch (i.e. existing behavior is unchanged).
      String? ref_label_column
      String? ref_batch_column

      # Genome for the ATAC cell-by-bin -> gene-activity conversion (multiome mode only):
      # "hg38" (default, human), "mm10", or "mm39" (mouse).
      String genome = "hg38"

      # Optional: also emit the per-cell maximum SCANVI posterior probability (the
      # confidence of the assigned label) as a `max_probability` obs column in every
      # output h5ad. Default false preserves existing outputs.
      Boolean output_max_probability = false

      # Optional pre-trained SCANVI model (.tar.gz of a saved model directory, no bundled adata).
      # When provided, the label-transfer task LOADS it and predicts instead of training SCVI/SCANVI,
      # and every query is predicted with the same model. Without it, each query trains its own
      # SCVI/SCANVI model on the reference plus that query. Every query also emits its model in
      # `scanvi_models_out` (same format), so an output feeds straight back here. Reuse is valid when
      # the model matches the incoming data/reference; novel-query re-mapping (scArches surgery) is
      # future work.
      File? scanvi_model

      # MultiomeLabelTransfer compute. Defaults give the GPU + a large box (training, or inference on
      # large data). Override (e.g. gpu_count=0 + small mem/cpu/disk) for a cheap CPU-only prediction
      # run, such as the pretrained Plumbing test.
      Int gpu_count = 2
      Int mem_size = 120
      Int nthreads = 32
      Int disk_size = 500

  }

  String pipeline_version = "3.0.0"

  # Docker image (same container for every task; only the label-transfer task gets GPUs attached)
  # Exposes run_gex_only_model for GEX-only mode (warp-tools/3rd-party-tools/scvi-scanvi).
  String docker = "us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:077d776c1957029b32ef79f12b6dac02e43d5673b70e915d1d690d95612a9831"

  # Multiome or GEX-only is decided here, once for every query, by whether ATAC files were given.
  Boolean multiome = length(atac_h5ads) > 0
  Int n_queries = length(gex_h5ads)
  Boolean inputs_paired = n_queries > 0 && length(input_ids) == n_queries && (!multiome || length(atac_h5ads) == n_queries)

  if (!inputs_paired) {
    call utils.ErrorWithMessage as ErrorInputsNotPaired {
        input:
          message = "input_ids (~{length(input_ids)}) and gex_h5ads (~{n_queries}) must be non-empty and the same length, and atac_h5ads (~{length(atac_h5ads)}) must be empty (GEX-only) or the same length (multiome). On Terra, set members with no atac_h5ad value shorten atac_h5ads, so run multiome and GEX-only samples as separate sets."
    }
  }

  # Step 1: CPU-only reference preparation, once for every query
  call PreprocessReference {
      input:
        ref_h5ad = ref_h5ad,
        ref_label_column = ref_label_column,
        ref_batch_column = ref_batch_column,
        docker = docker
  }

  # Steps 2 and 3 run once per query, independently of the other queries. Unpaired inputs give
  # no shards, so no query work starts before ErrorInputsNotPaired fails the run.
  scatter (i in range(if inputs_paired then n_queries else 0)) {
    if (multiome) {
      File atac_h5ad = atac_h5ads[i]
    }

    # Step 2: CPU-only preprocessing and filtering of this query's GEX (and ATAC)
    call PreprocessFilter {
        input:
          gex_h5ad = gex_h5ads[i],
          atac_h5ad = atac_h5ad,
          input_id = input_ids[i],
          genome = genome,
          docker = docker
    }

    # Step 3: SCVI/SCANVI model training + label transfer. Both tasks run the SAME container command
    # (label_transfer_from_preprocessed.py); they differ ONLY in the runtime GPU attributes, which
    # Cromwell cannot conditionally omit from a single task. So route by gpu_count: > 0 -> the GPU task,
    # 0 -> the CPU-only task (e.g. the pretrained-model prediction path). scvi-tools auto-detects the
    # accelerator, so the code path is identical either way.
    if (gpu_count > 0) {
      call MultiomeLabelTransfer {
          input:
            gex_h5ad = PreprocessFilter.preprocessed_gex_h5ad,
            atac_activity_h5ad = PreprocessFilter.preprocessed_atac_activity_h5ad,
            ref_h5ad = PreprocessReference.preprocessed_ref_h5ad,
            input_id = input_ids[i],
            max_epochs = max_epochs,
            batch_size = batch_size,
            seed = seed,
            output_max_probability = output_max_probability,
            scanvi_model = scanvi_model,
            gpu_count = gpu_count,
            mem_size = mem_size,
            nthreads = nthreads,
            disk_size = disk_size,
            docker = docker
      }
    }
    if (gpu_count == 0) {
      call MultiomeLabelTransferCpu {
          input:
            gex_h5ad = PreprocessFilter.preprocessed_gex_h5ad,
            atac_activity_h5ad = PreprocessFilter.preprocessed_atac_activity_h5ad,
            ref_h5ad = PreprocessReference.preprocessed_ref_h5ad,
            input_id = input_ids[i],
            max_epochs = max_epochs,
            batch_size = batch_size,
            seed = seed,
            output_max_probability = output_max_probability,
            scanvi_model = scanvi_model,
            mem_size = mem_size,
            nthreads = nthreads,
            disk_size = disk_size,
            docker = docker
      }
    }

    # Exactly one of the two tasks runs; select_first picks whichever produced each output.
    File scanvi_predictions_h5ad = select_first([MultiomeLabelTransfer.scanvi_predictions_h5ad, MultiomeLabelTransferCpu.scanvi_predictions_h5ad])
    File gex_annotated_h5ad = select_first([MultiomeLabelTransfer.gex_annotated_h5ad, MultiomeLabelTransferCpu.gex_annotated_h5ad])
    File scanvi_model_out = select_first([MultiomeLabelTransfer.scanvi_model_out, MultiomeLabelTransferCpu.scanvi_model_out])
    # The ATAC-annotated matrix exists in multiome mode only, where every query must produce one.
    if (multiome) {
      File atac_annotated_h5ad = select_first([MultiomeLabelTransfer.atac_annotated_h5ad, MultiomeLabelTransferCpu.atac_annotated_h5ad])
    }
  }

  output {
      # One entry per query, in input order. atac_annotated_h5ads is empty in GEX-only mode.
      Array[File] scanvi_predictions_h5ads = scanvi_predictions_h5ad
      Array[File] gex_annotated_h5ads = gex_annotated_h5ad
      Array[File] atac_annotated_h5ads = select_all(atac_annotated_h5ad)
      Array[File] scanvi_models_out = scanvi_model_out
      String pipeline_version_out = pipeline_version
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task 1: PreprocessReference (CPU-only, once per run)
#
# Adapts the annotated reference into the form the trainer expects, once for every query:
#   - Materializes counts from .raw for AIT-schema references (which have no .X)
#   - Maps the chosen label/batch obs columns to final_annotation/batch
#   - Tags the reference cells with modality rna_annotated
# ──────────────────────────────────────────────────────────────────────────────
task PreprocessReference {
    input {
        File ref_h5ad

        # Reference adaptation (AIT-schema support)
        String? ref_label_column
        String? ref_batch_column

        # Runtime attributes hardcoded in each task
        String docker
        Int disk_size = 1000
        Int mem_size = 120
        Int nthreads = 32
    }

    parameter_meta {
        ref_h5ad: "Annotated reference AnnData h5ad file. PBMC-style references carry cell type labels in obs['final_annotation'] and a batch in obs['batch']; AIT-schema references (uns['schema_version']+uns['hierarchy']) are auto-detected and adapted (see ref_label_column/ref_batch_column)."
        ref_label_column: "Reference obs column to use as the cell-type label. When unset, defaults to 'subclass' for AIT references and 'final_annotation' otherwise."
        ref_batch_column: "Reference obs column to use as the batch. When unset, defaults to 'donor_id' for AIT references and 'batch' otherwise."
        docker: "Docker image containing the scvi-scanvi runtime environment."
        disk_size: "Disk size in GB."
        mem_size: "Memory size in GB."
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import anndata as ad
import scanpy as sc

# Reference-adaptation controls (from WDL inputs; empty = use defaults).
ref_label_column = "~{default='' ref_label_column}".strip()
ref_batch_column = "~{default='' ref_batch_column}".strip()


def normalize_reference(ref, label_col_override, batch_col_override):
    """Adapt the annotated reference into the form the trainer expects: counts in
    .X, cell-type labels in obs['final_annotation'], and a batch in obs['batch'].

    Supports AIT-schema atlases (Allen Institute Taxonomy: no .X — counts live in
    .raw; labels in a class/subclass/cluster_id hierarchy; batch in donor_id) and
    PBMC-style references (already have .X + final_annotation + batch). For the
    latter with no overrides this is a validate-only no-op, so existing behavior
    (and truth) is unchanged.
    """
    is_ait = ("schema_version" in ref.uns) and ("hierarchy" in ref.uns)
    if is_ait:
        print(f"  AIT reference detected (schema_version={ref.uns.get('schema_version')}, "
              f"hierarchy levels={list(ref.uns['hierarchy'])}).")

    # AIT files store counts in .raw and have no .X; .raw.var loses the gene symbols,
    # so take the matrix from .raw but keep the (aligned) gene-symbol .var.
    if ref.X is None:
        if ref.raw is None:
            raise SystemExit("ERROR: reference has neither .X nor .raw; cannot obtain counts.")
        if ref.raw.n_vars != ref.n_vars:
            raise SystemExit(
                f"ERROR: reference .raw n_vars ({ref.raw.n_vars}) != .var n_vars ({ref.n_vars}); "
                "cannot align raw counts to gene symbols.")
        print("  Reference has no .X; materializing counts from .raw (with .var gene symbols).")
        counts = ref.raw.X
        if hasattr(counts, "tocsr"):
            counts = counts.tocsr()
        ref = ad.AnnData(X=counts, obs=ref.obs.copy(), var=ref.var.copy(), uns=dict(ref.uns))

    # Resolve label/batch columns: explicit override > AIT default > PBMC default.
    label_col = label_col_override or ("subclass" if is_ait else "final_annotation")
    batch_col = batch_col_override or ("donor_id" if is_ait else "batch")
    for col, role in [(label_col, "label"), (batch_col, "batch")]:
        if col not in ref.obs.columns:
            raise SystemExit(
                f"ERROR: reference {role} column '{col}' not found in obs "
                f"(available: {sorted(ref.obs.columns)}).")

    # Map into the trainer's expected columns. Skip when already correctly named so a
    # PBMC reference is left byte-for-byte unchanged.
    if label_col != "final_annotation":
        ref.obs["final_annotation"] = ref.obs[label_col].astype(str)
    if batch_col != "batch":
        ref.obs["batch"] = ref.obs[batch_col].astype(str)
    print(f"  Reference ready: label='{label_col}' -> final_annotation, "
          f"batch='{batch_col}' -> batch, shape {ref.shape}")
    return ref


print("Loading reference...")
ref = sc.read_h5ad("~{ref_h5ad}")
print(f"  Reference loaded: {ref.shape}")
# Adapt AIT-schema (or other) references into the trainer's expected form.
ref = normalize_reference(ref, ref_label_column, ref_batch_column)
ref.obs["modality"] = "rna_annotated"

print("Writing preprocessed reference...")
ref.write_h5ad("preprocessed_ref.h5ad")
print("Reference preprocessing complete.")
CODE
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        memory: "${mem_size} GiB"
        cpu: nthreads
        maxRetries: 1
    }

    output {
        File preprocessed_ref_h5ad = "preprocessed_ref.h5ad"
    }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task 2: PreprocessFilter (CPU-only, once per query)
#
# Handles one query's GEX (and, in multiome mode, ATAC) preprocessing and filtering
# before model training:
#   - Patches missing columns (star_IsCell, gex_barcodes)
#   - Filters GEX to STARsolo cell calls and min count thresholds
#   - Reindexes ATAC barcodes to match GEX
#   - Subsets to shared barcodes across GEX and ATAC
#   - Assigns batch labels and modality tags
#   - Converts ATAC cell-by-bin matrix to gene activity matrix (snapatac2)
#   - Outputs the preprocessed GEX (and ATAC activity) h5ad files ready for model training
# ──────────────────────────────────────────────────────────────────────────────
task PreprocessFilter {
    input {
        # This query's h5ad files; the workflow passes ATAC in multiome mode only
        File gex_h5ad
        File? atac_h5ad

        # Genome for the ATAC gene-activity conversion
        String genome = "hg38"

        # Runtime attributes hardcoded in each task
        String input_id
        String docker
        Int disk_size = 1000 # bigger disk before cell filtering
        Int mem_size = 120
        Int nthreads = 32
    }

    parameter_meta {
        gex_h5ad: "Gene expression AnnData h5ad file from Multiome/Optimus pipeline output."
        atac_h5ad: "Optional ATAC cell-by-bin AnnData h5ad file from Multiome/PeakCalling pipeline output, paired with gex_h5ad. If omitted, the query is preprocessed for GEX-only mode (training/annotation from the reference using GEX alone)."
        genome: "Genome for ATAC gene-activity conversion in multiome mode: 'hg38' (default), 'mm10', or 'mm39'."
        input_id: "Unique identifier prepended to all output filenames."
        docker: "Docker image containing the scvi-scanvi runtime environment."
        disk_size: "Disk size in GB."
        mem_size: "Memory size in GB."
    }

    command <<<
        set -euo pipefail

        # Symlink the gene annotation file (needed by snapatac2 make_gene_matrix);
        # only relevant when ATAC is present, but harmless to create unconditionally.
        ln -sf /usr/local/gencode.v41.basic.annotation.gff3.gz .

        # ── Run preprocessing in Python ───────────────────────────────────────
        python3 <<CODE
import scanpy as sc

gex_path  = "~{gex_h5ad}"
atac_path = "~{default='' atac_h5ad}"

# ATAC-genome control (from the WDL input; empty = use the default).
genome = "~{genome}".strip() or "hg38"

# The workflow passes ATAC in multiome mode only. Without it we run in GEX-only mode:
# train/annotate from the reference using GEX alone, never touching ATAC.
atac_present = bool(atac_path)
print(f"Mode: {'multiome (GEX + ATAC)' if atac_present else 'GEX-only (no ATAC)'}")

# ── 1. Load GEX (always required) ─────────────────────────────────────────
print("Loading GEX...")
gex = sc.read_h5ad(gex_path)
print(f"  GEX loaded: {gex.shape}")

# ── 2. Patch missing GEX columns ─────────────────────────────────────────
if "star_IsCell" not in gex.obs.columns:
    print("  Patching: adding star_IsCell = True (column was missing)")
    gex.obs["star_IsCell"] = True

# ── 3. Filter GEX to STARsolo cell calls ─────────────────────────────────
print("Filtering GEX to star_IsCell == True...")
gex = gex[gex.obs["star_IsCell"] == True]
print(f"  After star_IsCell filter: {gex.shape}")

# Filter genes and cells with fewer than 3 total counts
sc.pp.filter_genes(gex, min_counts=3)
sc.pp.filter_cells(gex, min_counts=3)
print(f"  After min_counts filter: {gex.shape}")

# Add batch label and preserve raw counts
gex.obs["batch"] = 1
gex.layers["counts"] = gex.X.copy()

input_id = "~{input_id}"

if atac_present:
    # ── Multiome path: load and integrate ATAC ───────────────────────────
    import snapatac2 as snap

    print("Loading ATAC...")
    atac = snap.read(atac_path, backed=None)
    print(f"  ATAC loaded: {atac.shape}")

    if "gex_barcodes" not in atac.obs.columns:
        print("  Patching: adding gex_barcodes from obs index (column was missing)")
        atac.obs["gex_barcodes"] = atac.obs.index

    # Reindex ATAC barcodes to match GEX
    print("Reindexing ATAC barcodes to gex_barcodes...")
    atac.obs = atac.obs.set_index("gex_barcodes")

    # Subset to shared barcodes
    shared = atac.obs.index.intersection(gex.obs.index)
    print(f"Shared barcodes between GEX and ATAC: {len(shared)}")
    atac_shared = atac[atac.obs.index.isin(shared)].copy()
    gex_shared  = gex[gex.obs.index.isin(shared)].copy()
    print(f"  GEX shared: {gex_shared.shape}")
    print(f"  ATAC shared: {atac_shared.shape}")

    # Assign batch labels
    atac_shared.obs["batch"] = "pd-multiome_sci_atac"
    gex_shared.obs["batch"]  = "pd-multiome_sci_gex"

    # Add placeholder annotations for query datasets
    gex_shared.obs["final_annotation"] = "Unknown"

    # Convert ATAC cell-by-bin → gene activity matrix (genome-specific)
    print(f"Converting ATAC cell-by-bin to gene activity matrix ({genome})...")
    _genomes = {
        "hg38": getattr(snap.genome, "hg38", None),
        "mm10": getattr(snap.genome, "mm10", None),
        "mm39": getattr(snap.genome, "mm39", None) or getattr(snap.genome, "GRCm39", None),
    }
    gene_anno = _genomes.get(genome)
    if gene_anno is None:
        raise SystemExit(f"ERROR: unsupported/unavailable genome '{genome}' (expected hg38|mm10|mm39)")
    atac_activity = snap.pp.make_gene_matrix(atac_shared, gene_anno=gene_anno)
    print(f"  Gene activity matrix: {atac_activity.shape}")
    atac_activity.obs["final_annotation"] = "Unknown"

    # Tag modalities
    gex_shared.obs["modality"]    = "rna_unannotated"
    atac_activity.obs["modality"] = "atac_unannotated"

    # Write preprocessed outputs (GEX + ATAC activity)
    print(f"Writing preprocessed outputs (input_id={input_id})...")
    gex_shared.write_h5ad(f"{input_id}_preprocessed_gex.h5ad")
    atac_activity.write_h5ad(f"{input_id}_preprocessed_atac_activity.h5ad")
else:
    # ── GEX-only path: no ATAC, no shared-barcode subset ─────────────────
    gex_shared = gex.copy()
    gex_shared.obs["batch"]            = "pd-multiome_sci_gex"
    gex_shared.obs["final_annotation"] = "Unknown"
    gex_shared.obs["modality"]         = "rna_unannotated"

    # Write preprocessed outputs (GEX only; no ATAC activity file)
    print(f"Writing preprocessed outputs (input_id={input_id})...")
    gex_shared.write_h5ad(f"{input_id}_preprocessed_gex.h5ad")

print("Preprocessing complete.")
CODE
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        memory: "${mem_size} GiB"
        cpu: nthreads
        maxRetries: 1
    }

    output {
        File preprocessed_gex_h5ad            = "~{input_id}_preprocessed_gex.h5ad"
        File? preprocessed_atac_activity_h5ad = "~{input_id}_preprocessed_atac_activity.h5ad"
    }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task 3: MultiomeLabelTransfer — SCVI/SCANVI model training + label transfer.
#
# Runs label_transfer_from_preprocessed.py in the container on one query's preprocessed h5ads
# (obtain the SCANVI model — load a supplied one or train — transfer labels, write annotated
# matrices + predictions + minimal model; no sample metadata is added). Two variants that are
# IDENTICAL except for the runtime GPU attributes (Cromwell cannot conditionally omit them from a
# single task):
#   - MultiomeLabelTransfer     requests a GPU (gpu_count > 0)
#   - MultiomeLabelTransferCpu  is CPU-only (gpu_count = 0, e.g. the pretrained-model prediction path)
# scvi-tools auto-detects the accelerator, so the command is identical for both. KEEP the two tasks'
# input / command / output blocks IN SYNC.
# ──────────────────────────────────────────────────────────────────────────────
task MultiomeLabelTransfer {
    input {
        # Preprocessed h5ad files: the query's from PreprocessFilter, the reference from PreprocessReference
        File gex_h5ad
        File? atac_activity_h5ad
        File ref_h5ad

        # Runtime attributes
        String input_id
        Int? max_epochs
        Int batch_size = 128
        Int? seed
        Boolean output_max_probability = false
        # Optional pre-trained SCANVI model (.tar.gz of a saved model dir); when provided, the task
        # loads it and predicts instead of training SCVI/SCANVI.
        File? scanvi_model
        String docker
        Int gpu_count = 2
        Int disk_size = 500
        Int mem_size = 120
        Int nthreads = 32
    }

    parameter_meta {
        gex_h5ad: "Preprocessed gene expression h5ad file (filtered, batch-labeled, modality-tagged)."
        atac_activity_h5ad: "Optional preprocessed ATAC gene activity h5ad file (converted from cell-by-bin, batch-labeled, modality-tagged). If omitted, the task runs in GEX-only mode and trains/annotates from GEX + reference only."
        ref_h5ad: "Preprocessed reference h5ad file with cell type labels and modality tag."
        input_id: "Unique identifier prepended to all output filenames."
        max_epochs: "Optional cap on SCVI/SCANVI training epochs, applied to both multiome and GEX-only modes. When unset, the container default (500) is used."
        seed: "Optional random seed for SCVI/SCANVI training. Unset (default) keeps training nondeterministic."
        batch_size: "SCVI/SCANVI minibatch size (SGD). Default 128 (reproduces prior behavior). Lower it to fit a high-cardinality reference on a small GPU (activation memory scales with batch_size x number of labels); raise it on a large-VRAM cloud GPU."
        output_max_probability: "When true, also write a `max_probability` obs column (the per-cell maximum SCANVI posterior probability, i.e. the assigned label's confidence) to every output h5ad."
        scanvi_model: "Optional .tar.gz of a saved SCANVI model directory (no bundled adata). When provided, the model is loaded and used to predict (no SCVI/SCANVI training). Valid when the model matches the incoming data/reference. The run still emits its model as scanvi_model_out."
        docker: "Docker image containing the scvi-scanvi runtime environment."
        gpu_count: "GPUs to request. Default 2 (training, or inference on large data)."
        disk_size: "Disk size in GB."
        mem_size: "Memory size in GB."
    }

    command <<<
        set -euo pipefail
        python3 /usr/local/label_transfer_from_preprocessed.py \
            --gex ~{gex_h5ad} \
            --ref ~{ref_h5ad} \
            --input-id ~{input_id} \
            --batch-size ~{batch_size} \
            ~{if defined(atac_activity_h5ad) then "--atac " + select_first([atac_activity_h5ad]) else ""} \
            ~{if defined(scanvi_model) then "--scanvi-model " + select_first([scanvi_model]) else ""} \
            ~{if defined(max_epochs) then "--max-epochs " + select_first([max_epochs]) else ""} \
            ~{"--seed " + seed} \
            ~{true="--output-max-probability" false="" output_max_probability}
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        memory: "${mem_size} GiB"
        cpu: nthreads
        # camelCase GPU keys = portable across Cromwell/Terra/GCP (per WARP conventions).
        gpuType: "nvidia-tesla-t4"           # known to work with Terra
        gpuCount: gpu_count
        nvidiaDriverVersion: "535.104.05"    # compatible with CUDA 12.x and T4 GPUs
        maxRetries: 1
    }

    output {
        File scanvi_predictions_h5ad = "~{input_id}_SCANVI_predictions.h5ad"
        File? atac_annotated_h5ad    = "~{input_id}_atac_annotated_matrix.h5ad"
        File gex_annotated_h5ad      = "~{input_id}_gex_annotated_matrix.h5ad"
        File scanvi_model_out        = "~{input_id}_scanvi_model.tar.gz"
    }
}


# CPU-only variant of MultiomeLabelTransfer: identical inputs/command/outputs, but NO GPU runtime
# attributes and no gpu_count input. Selected by the workflow when gpu_count = 0 (e.g. the
# pretrained-model prediction path). Keep in sync with MultiomeLabelTransfer above.
task MultiomeLabelTransferCpu {
    input {
        # Preprocessed h5ad files: the query's from PreprocessFilter, the reference from PreprocessReference
        File gex_h5ad
        File? atac_activity_h5ad
        File ref_h5ad

        # Runtime attributes
        String input_id
        Int? max_epochs
        Int batch_size = 128
        Int? seed
        Boolean output_max_probability = false
        File? scanvi_model
        String docker
        Int disk_size = 500
        Int mem_size = 120
        Int nthreads = 32
    }

    parameter_meta {
        gex_h5ad: "Preprocessed gene expression h5ad file (filtered, batch-labeled, modality-tagged)."
        atac_activity_h5ad: "Optional preprocessed ATAC gene activity h5ad file. If omitted, GEX-only mode."
        ref_h5ad: "Preprocessed reference h5ad file with cell type labels and modality tag."
        input_id: "Unique identifier prepended to all output filenames."
        max_epochs: "Optional cap on SCVI/SCANVI training epochs. When unset, the container default (500) is used."
        batch_size: "SCVI/SCANVI minibatch size (SGD). Default 128."
        seed: "Optional random seed for SCVI/SCANVI training. Unset (default) keeps training nondeterministic."
        output_max_probability: "When true, also write a `max_probability` obs column to every output h5ad."
        scanvi_model: "Optional .tar.gz of a saved SCANVI model directory (no bundled adata). When provided, the model is loaded and used to predict (no SCVI/SCANVI training)."
        docker: "Docker image containing the scvi-scanvi runtime environment."
        disk_size: "Disk size in GB."
        mem_size: "Memory size in GB."
    }

    command <<<
        set -euo pipefail
        python3 /usr/local/label_transfer_from_preprocessed.py \
            --gex ~{gex_h5ad} \
            --ref ~{ref_h5ad} \
            --input-id ~{input_id} \
            --batch-size ~{batch_size} \
            ~{if defined(atac_activity_h5ad) then "--atac " + select_first([atac_activity_h5ad]) else ""} \
            ~{if defined(scanvi_model) then "--scanvi-model " + select_first([scanvi_model]) else ""} \
            ~{if defined(max_epochs) then "--max-epochs " + select_first([max_epochs]) else ""} \
            ~{"--seed " + seed} \
            ~{true="--output-max-probability" false="" output_max_probability}
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        memory: "${mem_size} GiB"
        cpu: nthreads
        maxRetries: 1
    }

    output {
        File scanvi_predictions_h5ad = "~{input_id}_SCANVI_predictions.h5ad"
        File? atac_annotated_h5ad    = "~{input_id}_atac_annotated_matrix.h5ad"
        File gex_annotated_h5ad      = "~{input_id}_gex_annotated_matrix.h5ad"
        File scanvi_model_out        = "~{input_id}_scanvi_model.tar.gz"
    }
}
