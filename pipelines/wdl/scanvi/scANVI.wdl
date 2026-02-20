version 1.0

import "../../../tasks/wdl/Utilities.wdl" as utils

workflow ScviScanvi {

  meta {
    description: "Pipeline for cell type label transfer on Multiome data using SCVI and SCANVI models. Integrates single-cell RNA (GEX) and ATAC data with an annotated reference to transfer cell type labels via semi-supervised deep generative models."
    allowNestedInputs: true
  }

  input {
      # Required h5ad inputs
      File gex_h5ad
      File atac_h5ad
      File ref_h5ad

      # Runtime attributes
      String cloud_provider
      Int disk_size = 500
      Int mem_size = 64
      Int nthreads = 8

      # GPU configuration
      String gpu_type = "nvidia-tesla-t4"
      Int gpu_count = 1
  }

  String pipeline_version = "1.0.0"

  # Determine docker prefix based on cloud provider
  String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
  String acr_docker_prefix = "dsppipelinedev.azurecr.io/"
  String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

  # Docker image
  String scvi_scanvi_docker = "scvi-scanvi:1.0.0-1.2-1756234975"

  # Make sure either 'gcp' or 'azure' is supplied as cloud_provider input. If not, raise an error
  if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
      call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
        input:
          message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
      }
  }

  call MultiomeLabelTransfer {
      input:
        gex_h5ad = gex_h5ad,
        atac_h5ad = atac_h5ad,
        ref_h5ad = ref_h5ad,
        docker_path = docker_prefix + scvi_scanvi_docker,
        disk_size = disk_size,
        mem_size = mem_size,
        nthreads = nthreads,
        gpu_type = gpu_type,
        gpu_count = gpu_count
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
        File gex_h5ad
        File atac_h5ad
        File ref_h5ad

        # Runtime attributes
        String docker_path
        Int disk_size = 500
        Int mem_size = 64
        Int nthreads = 8
        String gpu_type = "nvidia-tesla-t4"
        Int gpu_count = 1
    }

    parameter_meta {
        gex_h5ad: "Gene expression AnnData h5ad file from Multiome/Optimus pipeline output."
        atac_h5ad: "ATAC cell-by-bin AnnData h5ad file from Multiome/PeakCalling pipeline output."
        ref_h5ad: "Annotated reference AnnData h5ad file with cell type labels in obs['final_annotation']."
        docker_path: "Docker image path containing the scvi-scanvi runtime environment."
        disk_size: "Disk size in GB."
        mem_size: "Memory size in GB."
        gpu_type: "GPU type for accelerated model training."
        gpu_count: "Number of GPUs to use."
    }

    command <<<
        set -euo pipefail

        python3 /usr/local/multiome_label_transfer.py \
            --gex-file ~{gex_h5ad} \
            --atac-file ~{atac_h5ad} \
            --ref-file ~{ref_h5ad}
    >>>

    runtime {
        docker: docker_path
        disks: "local-disk ${disk_size} SSD"
        memory: "${mem_size} GiB"
        cpu: nthreads
        gpuType: gpu_type
        gpuCount: gpu_count
        zones: ["us-central1-c"]
    }

    output {
        File scanvi_predictions_h5ad = "SCANVI_predictions.h5ad"
        File atac_annotated_h5ad = "atac_annotated_matrix.h5ad"
        File gex_annotated_h5ad = "gex_annotated_matrix.h5ad"
    }
}
