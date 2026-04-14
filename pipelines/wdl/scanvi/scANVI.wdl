version 1.0

workflow scANVI {

  meta {
    description: "Pipeline for cell type label transfer on Multiome data using SCVI and SCANVI models. Integrates single-cell RNA (GEX) and ATAC data with an annotated reference to transfer cell type labels via semi-supervised deep generative models."
    allowNestedInputs: true
  }

  input {
      # GCS bucket path containing input h5ad files (e.g., gs://bucket/path/to/inputs)
      String? input_bucket

      # Optional direct h5ad file inputs (override bucket files if provided)
      File? gex_h5ad
      File? atac_h5ad
      File? ref_h5ad

      # Expected filenames in the input bucket
      String gex_filename = "gex.h5ad"
      String atac_filename = "atac.h5ad"
      String ref_filename = "ref.h5ad"

  }

  String pipeline_version = "1.0.1"

  # Docker image (same container for both tasks; only Task 2 gets GPUs attached)
  String docker = "us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:81fe915a045bd2929a1c457f4a0061055c6ea42fa3f88e9352b618e4a6e47b58"
  # Step 1: CPU-only preprocessing and filtering of all three h5ad inputs
  call PreprocessFilter {
      input:
        input_bucket = input_bucket,
        gex_h5ad = gex_h5ad,
        atac_h5ad = atac_h5ad,
        ref_h5ad = ref_h5ad,
        gex_filename = gex_filename,
        atac_filename = atac_filename,
        ref_filename = ref_filename,
        docker = docker
  }

  # Step 2: GPU-accelerated SCVI/SCANVI model training and label transfer
  call MultiomeLabelTransfer {
      input:
        gex_h5ad = PreprocessFilter.preprocessed_gex_h5ad,
        atac_activity_h5ad = PreprocessFilter.preprocessed_atac_activity_h5ad,
        ref_h5ad = PreprocessFilter.preprocessed_ref_h5ad,
        docker = docker
  }

  output {
      File scanvi_predictions_h5ad = MultiomeLabelTransfer.scanvi_predictions_h5ad
      File atac_annotated_h5ad = MultiomeLabelTransfer.atac_annotated_h5ad
      File gex_annotated_h5ad = MultiomeLabelTransfer.gex_annotated_h5ad
      String pipeline_version_out = pipeline_version
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task 1: PreprocessFilter (CPU-only)
#
# Handles all h5ad preprocessing and filtering before model training:
#   - Patches missing columns (star_IsCell, gex_barcodes)
#   - Filters GEX to STARsolo cell calls and min count thresholds
#   - Reindexes ATAC barcodes to match GEX
#   - Subsets to shared barcodes across GEX and ATAC
#   - Assigns batch labels and modality tags
#   - Converts ATAC cell-by-bin matrix to gene activity matrix (snapatac2)
#   - Outputs three preprocessed h5ad files ready for model training
# ──────────────────────────────────────────────────────────────────────────────
task PreprocessFilter {
    input {
        # GCS bucket path containing input h5ad files
        String? input_bucket

        # Optional direct h5ad file inputs
        File? gex_h5ad
        File? atac_h5ad
        File? ref_h5ad

        # Expected filenames in the input bucket
        String gex_filename = "gex.h5ad"
        String atac_filename = "atac.h5ad"
        String ref_filename = "ref.h5ad"

        # Runtime attributes hardcoded in each task
        String docker
        Int disk_size = 1000 # bigger disk before cell filtering
        Int mem_size = 120
        Int nthreads = 32
    }

    parameter_meta {
        input_bucket: "GCS bucket path containing input h5ad files (e.g., gs://bucket/path/to/inputs)."
        gex_h5ad: "Gene expression AnnData h5ad file from Multiome/Optimus pipeline output."
        atac_h5ad: "ATAC cell-by-bin AnnData h5ad file from Multiome/PeakCalling pipeline output."
        ref_h5ad: "Annotated reference AnnData h5ad file with cell type labels in obs['final_annotation']."
        gex_filename: "Expected GEX h5ad filename in the input bucket."
        atac_filename: "Expected ATAC h5ad filename in the input bucket."
        ref_filename: "Expected reference h5ad filename in the input bucket."
        docker: "Docker image containing the scvi-scanvi runtime environment."
        disk_size: "Disk size in GB."
        mem_size: "Memory size in GB."
    }

    command <<<
        set -euo pipefail

        # ── Resolve input file paths ──────────────────────────────────────────
        if [ -n "~{default='' gex_h5ad}" ]; then
            # Direct File inputs: Cromwell already localized them
            GEX_FILE="~{default='' gex_h5ad}"
            ATAC_FILE="~{default='' atac_h5ad}"
            REF_FILE="~{default='' ref_h5ad}"
        else
            # Bucket mode: construct GCS paths and download
            BUCKET="~{default='' input_bucket}"
            BUCKET="${BUCKET%/}"
            GEX_FILE="${BUCKET}/~{gex_filename}"
            ATAC_FILE="${BUCKET}/~{atac_filename}"
            REF_FILE="${BUCKET}/~{ref_filename}"

            echo "Downloading inputs from bucket..."
            gsutil cp "$GEX_FILE" gex_input.h5ad
            gsutil cp "$ATAC_FILE" atac_input.h5ad
            gsutil cp "$REF_FILE" ref_input.h5ad
            GEX_FILE="gex_input.h5ad"
            ATAC_FILE="atac_input.h5ad"
            REF_FILE="ref_input.h5ad"
        fi

        export GEX_FILE ATAC_FILE REF_FILE

        # Symlink the gene annotation file (needed by snapatac2 make_gene_matrix)
        ln -sf /usr/local/gencode.v41.basic.annotation.gff3.gz .

        # ── Run preprocessing in Python ───────────────────────────────────────
        python3 <<'PYEOF'
import os
import anndata as ad
import scanpy as sc
import snapatac2 as snap
import numpy as np

gex_path  = os.environ["GEX_FILE"]
atac_path = os.environ["ATAC_FILE"]
ref_path  = os.environ["REF_FILE"]

# ── 1. Load datasets ─────────────────────────────────────────────────────
print("Loading GEX...")
gex = sc.read_h5ad(gex_path)
print(f"  GEX loaded: {gex.shape}")

print("Loading ATAC...")
atac = snap.read(atac_path, backed=None)
print(f"  ATAC loaded: {atac.shape}")

print("Loading reference...")
ref = sc.read_h5ad(ref_path)
print(f"  Reference loaded: {ref.shape}")

# ── 2. Patch missing columns ─────────────────────────────────────────────
if "star_IsCell" not in gex.obs.columns:
    print("  Patching: adding star_IsCell = True (column was missing)")
    gex.obs["star_IsCell"] = True

if "gex_barcodes" not in atac.obs.columns:
    print("  Patching: adding gex_barcodes from obs index (column was missing)")
    atac.obs["gex_barcodes"] = atac.obs.index

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

# ── 4. Reindex ATAC barcodes to match GEX ────────────────────────────────
print("Reindexing ATAC barcodes to gex_barcodes...")
atac.obs = atac.obs.set_index("gex_barcodes")

# ── 5. Subset to shared barcodes ─────────────────────────────────────────
shared = atac.obs.index.intersection(gex.obs.index)
print(f"Shared barcodes between GEX and ATAC: {len(shared)}")
atac_shared = atac[atac.obs.index.isin(shared)].copy()
gex_shared  = gex[gex.obs.index.isin(shared)].copy()
print(f"  GEX shared: {gex_shared.shape}")
print(f"  ATAC shared: {atac_shared.shape}")

# ── 6. Assign batch labels ───────────────────────────────────────────────
atac_shared.obs["batch"] = "pd-multiome_sci_atac"
gex_shared.obs["batch"]  = "pd-multiome_sci_gex"

# ── 7. Add placeholder annotations for query datasets ────────────────────
gex_shared.obs["final_annotation"] = "Unknown"

# ── 8. Convert ATAC cell-by-bin → gene activity matrix ───────────────────
print("Converting ATAC cell-by-bin to gene activity matrix (hg38)...")
atac_activity = snap.pp.make_gene_matrix(atac_shared, gene_anno=snap.genome.hg38)
print(f"  Gene activity matrix: {atac_activity.shape}")
atac_activity.obs["final_annotation"] = "Unknown"

# ── 9. Tag modalities ────────────────────────────────────────────────────
gex_shared.obs["modality"]    = "rna_unannotated"
atac_activity.obs["modality"] = "atac_unannotated"
ref.obs["modality"]           = "rna_annotated"

# ── 10. Write preprocessed outputs ───────────────────────────────────────
print("Writing preprocessed outputs...")
gex_shared.write_h5ad("preprocessed_gex.h5ad")
atac_activity.write_h5ad("preprocessed_atac_activity.h5ad")
ref.write_h5ad("preprocessed_ref.h5ad")
print("Preprocessing complete.")
PYEOF
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
        File preprocessed_gex_h5ad           = "preprocessed_gex.h5ad"
        File preprocessed_atac_activity_h5ad = "preprocessed_atac_activity.h5ad"
        File preprocessed_ref_h5ad           = "preprocessed_ref.h5ad"
    }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task 2: MultiomeLabelTransfer (GPU)
#
# Receives preprocessed h5ad files and runs model training + label transfer:
#   - SCVI unsupervised latent space learning
#   - SCANVI semi-supervised label transfer
#   - Outputs annotated GEX, ATAC, and SCANVI prediction h5ad files
# ──────────────────────────────────────────────────────────────────────────────
task MultiomeLabelTransfer {
    input {
        # Preprocessed h5ad files from PreprocessFilter
        File gex_h5ad
        File atac_activity_h5ad
        File ref_h5ad

        # Runtime attributes
        String docker
        Int disk_size = 500
        Int mem_size = 120
        Int nthreads = 32
    }

    parameter_meta {
        gex_h5ad: "Preprocessed gene expression h5ad file (filtered, batch-labeled, modality-tagged)."
        atac_activity_h5ad: "Preprocessed ATAC gene activity h5ad file (converted from cell-by-bin, batch-labeled, modality-tagged)."
        ref_h5ad: "Preprocessed reference h5ad file with cell type labels and modality tag."
        docker: "Docker image containing the scvi-scanvi runtime environment."
        disk_size: "Disk size in GB."
        mem_size: "Memory size in GB."
    }

    command <<<
        set -euo pipefail

        python3 <<'PYEOF'
import sys
import time
import pandas as pd
import scanpy as sc
import anndata as ad

# ── Import ONLY model-training and label-transfer functions ───────────────
# We import specific functions from the container script so that main()
# (which re-runs all preprocessing) is NEVER called.  The three functions
# imported here operate entirely on already-preprocessed AnnData objects.
sys.path.insert(0, "/usr/local")
from multiome_label_transfer import run_multi_model, transfer_labels, finalize_output

# ── 1. Load preprocessed h5ad files ──────────────────────────────────────
# These files were fully preprocessed by the PreprocessFilter task:
#   - GEX: filtered to STARsolo cell calls, min-count filtered, batch-
#     labeled, modality-tagged, subset to shared barcodes
#   - ATAC activity: gene activity matrix (converted from cell-by-bin by
#     snapatac2), batch-labeled, modality-tagged, subset to shared barcodes
#   - Reference: annotated with final_annotation and modality tag
# NO additional filtering, reindexing, or gene-activity conversion needed.
print("Loading preprocessed inputs (no re-preprocessing)...")
gex           = sc.read_h5ad("~{gex_h5ad}")
atac_activity = sc.read_h5ad("~{atac_activity_h5ad}")
ref           = sc.read_h5ad("~{ref_h5ad}")
print(f"  GEX:           {gex.shape}")
print(f"  ATAC activity: {atac_activity.shape}")
print(f"  Reference:     {ref.shape}")

timing_summary = {}

# ── 2. Train SCVI and SCANVI models (GPU-accelerated) ────────────────────
# run_multi_model() performs the following steps on the preprocessed data:
#   a. Concatenates GEX, ATAC activity, and reference into a single AnnData
#      (ad.concat with join='inner', index_unique='_')
#   b. Filters genes expressed in fewer than 5 cells
#   c. Selects 5000 highly variable genes (Seurat v3, batch-aware)
#   d. Trains SCVI (unsupervised VAE): 2 layers, 30 latent dims,
#      negative-binomial likelihood, gene-batch dispersion, up to 500 epochs
#   e. Trains SCANVI (semi-supervised): propagates reference cell-type labels
#      to unlabeled query cells, up to 500 epochs, 100 samples per label
#   f. Returns: concatenated AnnData (data), SCVI model (vae), SCANVI (lvae)
print("Training SCVI and SCANVI models...")
start = time.time()
data, vae, lvae = run_multi_model(gex, atac_activity, ref)
timing_summary['Model Training'] = time.time() - start
print(f"  Model training complete in {timing_summary['Model Training']:.1f}s")

# ── 3. Transfer labels using the trained SCANVI model ────────────────────
# transfer_labels() performs the following steps:
#   a. Predicts cell-type labels (C_scANVI) for every cell using SCANVI
#   b. Extracts the SCANVI latent representation (X_scANVI)
#   c. Computes a neighborhood graph and UMAP from the latent space
#   d. Writes intermediate labeled AnnData (adata_scanvi_labels.h5ad)
#   e. Returns a UMAP plot and the updated AnnData object
print("Transferring labels with SCANVI...")
start = time.time()
plot, data = transfer_labels(data, lvae)
timing_summary['Label Transfer'] = time.time() - start
print(f"  Label transfer complete in {timing_summary['Label Transfer']:.1f}s")

# ── 4. Propagate predicted labels back to original matrices ──────────────
# ad.concat (called inside run_multi_model) appended modality-key suffixes
# to obs_names (e.g. "barcode_rna_unannotated", "barcode_atac_unannotated").
# We use those suffixed names to look up each cell's SCANVI prediction
# (C_scANVI) in the concatenated object and copy it into the original
# preprocessed GEX and ATAC AnnData objects.
print("Propagating predicted labels to GEX and ATAC matrices...")
gex.obs['final_annotation'] = data.obs.loc[
    gex.obs_names + '_rna_unannotated']['C_scANVI'].to_numpy()
atac_activity.obs['final_annotation'] = data.obs.loc[
    atac_activity.obs_names + '_atac_unannotated']['C_scANVI'].to_numpy()

# ── 5. Write annotated GEX and ATAC matrices ────────────────────────────
# These outputs mirror what main() in the container script produces,
# but are built from the already-preprocessed inputs.
print("Writing annotated matrices...")
gex.write("gex_annotated_matrix.h5ad")
atac_activity.write("atac_annotated_matrix.h5ad")
print(f"  gex_annotated_matrix.h5ad:  {gex.shape}")
print(f"  atac_annotated_matrix.h5ad: {atac_activity.shape}")

# ── 6. Finalize SCANVI predictions ──────────────────────────────────────
# finalize_output() adds downstream-required metadata and reformats:
#   a. Adds placeholder metadata (biosample_id, donor_id, species,
#      disease, organ, library_preparation_protocol, sex)
#   b. Renames 'final_annotation' → 'celltype'
#   c. Copies the count matrix into the .raw layer for SCP ingest
print("Finalizing SCANVI predictions...")
final_data = finalize_output(data)
final_data.write("SCANVI_predictions.h5ad")
print(f"  SCANVI_predictions.h5ad:    {final_data.shape}")

# ── Timing summary ───────────────────────────────────────────────────────
timing_df = pd.DataFrame(
    [(step, f"{elapsed:.2f}s") for step, elapsed in timing_summary.items()],
    columns=["Step", "Time"]
)
print("\nTiming Summary:")
print(timing_df.to_string(index=False))
print("\nMultiomeLabelTransfer complete.")
PYEOF
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        memory: "${mem_size} GiB"
        cpu: nthreads
        gpuType: "nvidia-tesla-t4"
        gpuCount: 2
        nvidiaDriverVersion: "525.147.05"
        zones: "us-central1-a us-central1-c"
        maxRetries: 1
    }

    output {
        File scanvi_predictions_h5ad = "SCANVI_predictions.h5ad"
        File atac_annotated_h5ad = "atac_annotated_matrix.h5ad"
        File gex_annotated_h5ad = "gex_annotated_matrix.h5ad"
    }
}
