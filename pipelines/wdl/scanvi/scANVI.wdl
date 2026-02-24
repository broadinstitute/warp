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
  String scvi_scanvi_docker = "scvi-scanvi:1.0.0-1.2-1760025671"

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

        # Download files from input bucket if provided
        if [ -n "~{default='' input_bucket}" ]; then
            python3 -c "
import sys
sys.path.insert(0, '/usr/local')
from gcs_utils import pull_all_files
bucket_path = '~{default='' input_bucket}'.rstrip('/')
files = [
    bucket_path + '/~{gex_filename}',
    bucket_path + '/~{atac_filename}',
    bucket_path + '/~{ref_filename}'
]
pull_all_files(files)
"
        fi

        # Use direct file inputs if provided, otherwise use bucket-downloaded files
        if [ -n "~{default='' gex_h5ad}" ]; then
            GEX_FILE="~{default='' gex_h5ad}"
        else
            GEX_FILE="~{gex_filename}"
        fi

        if [ -n "~{default='' atac_h5ad}" ]; then
            ATAC_FILE="~{default='' atac_h5ad}"
        else
            ATAC_FILE="~{atac_filename}"
        fi

        if [ -n "~{default='' ref_h5ad}" ]; then
            REF_FILE="~{default='' ref_h5ad}"
        else
            REF_FILE="~{ref_filename}"
        fi

        python3 /usr/local/multiome_label_transfer.py \
            --gex-file "$GEX_FILE" \
            --atac-file "$ATAC_FILE" \
            --ref-file "$REF_FILE"
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
