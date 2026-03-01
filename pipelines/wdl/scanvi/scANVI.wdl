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

      # Runtime attributes
      String cloud_provider = "gcp"
      Int disk_size = 500
      Int mem_size = 120
      Int nthreads = 32

      # GPU configuration
      String gpu_type = "nvidia-tesla-t4"
      Int gpu_count = 2
      String nvidiaDriverVersion = "535.104.05"
  }

  String pipeline_version = "1.0.0"

  # Docker image
  String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
  String docker_prefix = gcr_docker_prefix
  String scvi_scanvi_docker = "scvi-scanvi:rc_3220_scanvi"

  call MultiomeLabelTransfer {
      input:
        input_bucket = input_bucket,
        gex_h5ad = gex_h5ad,
        atac_h5ad = atac_h5ad,
        ref_h5ad = ref_h5ad,
        gex_filename = gex_filename,
        atac_filename = atac_filename,
        ref_filename = ref_filename,
        docker_path = docker_prefix + scvi_scanvi_docker,
        disk_size = disk_size,
        mem_size = mem_size,
        nthreads = nthreads,
        gpu_type = gpu_type,
        gpu_count = gpu_count,
        nvidiaDriverVersion = nvidiaDriverVersion
  }

  output {
      File scanvi_predictions_h5ad = MultiomeLabelTransfer.scanvi_predictions_h5ad
      File atac_annotated_h5ad = MultiomeLabelTransfer.atac_annotated_h5ad
      File gex_annotated_h5ad = MultiomeLabelTransfer.gex_annotated_h5ad
      String pipeline_version_out = pipeline_version
  }
}

task MultiomeLabelTransfer {
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

        # Runtime attributes
        String docker_path
        Int disk_size = 500
        Int mem_size = 120
        Int nthreads = 32
        String gpu_type = "nvidia-tesla-t4"
        Int gpu_count = 2
        String nvidiaDriverVersion = "535.104.05"
    }

    parameter_meta {
        input_bucket: "GCS bucket path containing input h5ad files (e.g., gs://bucket/path/to/inputs)."
        gex_h5ad: "Gene expression AnnData h5ad file from Multiome/Optimus pipeline output."
        atac_h5ad: "ATAC cell-by-bin AnnData h5ad file from Multiome/PeakCalling pipeline output."
        ref_h5ad: "Annotated reference AnnData h5ad file with cell type labels in obs['final_annotation']."
        gex_filename: "Expected GEX h5ad filename in the input bucket."
        atac_filename: "Expected ATAC h5ad filename in the input bucket."
        ref_filename: "Expected reference h5ad filename in the input bucket."
        docker_path: "Docker image path containing the scvi-scanvi runtime environment."
        disk_size: "Disk size in GB."
        mem_size: "Memory size in GB."
        gpu_type: "GPU type for accelerated model training."
        gpu_count: "Number of GPUs to use."
        nvidiaDriverVersion: "NVIDIA driver version for GPU support."
    }

    command <<<
        set -euo pipefail

        # Build file paths and determine localize mode
        LOCALIZE_FLAG=""

        if [ -n "~{default='' gex_h5ad}" ]; then
            # Direct File inputs: Cromwell already localized them
            GEX_FILE="~{default='' gex_h5ad}"
            ATAC_FILE="~{default='' atac_h5ad}"
            REF_FILE="~{default='' ref_h5ad}"
        else
            # Bucket mode: construct GCS paths and let the script download them
            BUCKET="~{default='' input_bucket}"
            BUCKET="${BUCKET%/}"
            GEX_FILE="${BUCKET}/~{gex_filename}"
            ATAC_FILE="${BUCKET}/~{atac_filename}"
            REF_FILE="${BUCKET}/~{ref_filename}"
            LOCALIZE_FLAG="--localize"
        fi

        # Ensure GEX h5ad has 'star_IsCell' column and ATAC h5ad has 'gex_barcodes'
        # column (both required by the script). If missing, add sensible defaults.
        python3 -c "
import anndata as ad
import os

# Patch GEX: add star_IsCell if missing (makes the cell filter a no-op)
gex_path = '$GEX_FILE'
if not gex_path.startswith('gs://') and os.path.exists(gex_path):
    gex = ad.read_h5ad(gex_path)
    if 'star_IsCell' not in gex.obs.columns:
        print('Adding missing star_IsCell column (all True) to GEX h5ad')
        gex.obs['star_IsCell'] = True
        gex.write_h5ad(gex_path)
        print('Patched GEX h5ad saved')
    else:
        print('star_IsCell column already present')

# Patch ATAC: add gex_barcodes if missing (uses existing obs index as barcodes)
atac_path = '$ATAC_FILE'
if not atac_path.startswith('gs://') and os.path.exists(atac_path):
    atac = ad.read_h5ad(atac_path)
    if 'gex_barcodes' not in atac.obs.columns:
        print('Adding missing gex_barcodes column (copy of obs index) to ATAC h5ad')
        atac.obs['gex_barcodes'] = atac.obs.index
        atac.write_h5ad(atac_path)
        print('Patched ATAC h5ad saved')
    else:
        print('gex_barcodes column already present')
"

        # Symlink the gene annotation file to the working directory
        # (the script uses a relative path to find it)
        ln -sf /usr/local/gencode.v41.basic.annotation.gff3.gz .

        python3 /usr/local/multiome_label_transfer.py \
            --gex-file "$GEX_FILE" \
            --atac-file "$ATAC_FILE" \
            --ref-file "$REF_FILE" \
            $LOCALIZE_FLAG
    >>>

    runtime {
        docker: docker_path
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        memory: "${mem_size} GiB"
        cpu: nthreads
        gpuType: gpu_type
        gpuCount: gpu_count
        nvidiaDriverVersion: nvidiaDriverVersion
        zones: ["us-central1-a", "us-central1-c"]
        maxRetries: 1
    }

    output {
        File scanvi_predictions_h5ad = "SCANVI_predictions.h5ad"
        File atac_annotated_h5ad = "atac_annotated_matrix.h5ad"
        File gex_annotated_h5ad = "gex_annotated_matrix.h5ad"
    }
}
