version 1.0

import "../../../tasks/wdl/Utilities.wdl" as utils

workflow scANVI {

  meta {
    description: "Pipeline for cell type label transfer using SCVI and SCANVI models. Labels one or more samples, each independently and with its own outputs, either with a supplied SCANVI model or by training on an annotated reference plus the sample. Each sample's single-cell RNA (GEX) and, optionally, ATAC data are labelled; the annotated matrices are written ready for Single Cell Portal. All samples of a run go through one machine, one after another."
    allowNestedInputs: true
  }

  input {
      # Samples, paired by position: input_ids[i] prefixes every output of gex_h5ads[i]
      # (and atac_h5ads[i]) and is recorded in each output's uns. A single sample is a
      # one-element array.
      Array[String] input_ids
      Array[File] gex_h5ads

      # ATAC cell-by-bin h5ads for multiome mode: one per gex_h5ads entry, in the same order.
      # Empty (default) runs every sample in GEX-only mode.
      Array[File] atac_h5ads = []

      # Optional pre-trained SCANVI model (.tar.gz of a saved model directory, no bundled adata),
      # used for every sample. A model this pipeline trained loads as is; a model trained on the
      # reference alone (warp-scanvi-training) gets each sample mapped onto it (scArches).
      # Without a model, each sample trains its own SCVI/SCANVI model on ref_h5ad plus that sample.
      File? scanvi_model

      # Annotated reference h5ad. Needed only to train, i.e. when no scanvi_model is given.
      File? ref_h5ad

      # Optional cap on SCVI/SCANVI training epochs. When unset, the container default (500) is used.
      Int? max_epochs

      # Epochs when mapping a sample onto a reference-only model. When unset, the container
      # default (100) is used.
      Int? query_max_epochs

      # SCVI/SCANVI minibatch size (SGD). Default 128 (scvi-tools' own default). Lower it (e.g. 8)
      # to fit a high-cardinality reference on a small GPU; raise it on a large-VRAM cloud GPU.
      Int batch_size = 128

      # Optional random seed for SCVI/SCANVI training. Unset (default) keeps training
      # nondeterministic; set it to make a run reproducible (tests pin one).
      Int? seed

      # Reference adaptation when training. For an AIT-schema atlas these select which obs columns
      # become the cell-type label and the batch. When unset, AIT references default to
      # subclass/donor_id and PBMC-style references to final_annotation/batch.
      String? ref_label_column
      String? ref_batch_column

      # Genome for the ATAC cell-by-bin -> gene-activity conversion (multiome mode only):
      # "hg38" (default, human), "mm10", or "mm39" (mouse).
      String genome = "hg38"

      # Optional: also emit the per-cell maximum SCANVI posterior probability (the confidence of
      # the assigned label) as a `max_probability` obs column in every output h5ad.
      Boolean output_max_probability = false

      # Label-transfer machine. Every sample runs on it in turn, so memory only has to fit the
      # largest sample. The disk is sized from the inputs unless set. For a supplied-model run on
      # CPU, set gpu_count = 0.
      Int gpu_count = 2
      Int mem_size = 104
      Int nthreads = 16
      Int? disk_size

  }

  String pipeline_version = "3.0.0"

  # Docker image (same container for every task; only the label-transfer task gets GPUs attached)
  String docker = "us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:5764b99116daedef6c73a83483028ac78e318b1ec33de8824a9861b4d3201b52"

  # Multiome or GEX-only is decided here, once for every sample, by whether ATAC files were given.
  Boolean multiome = length(atac_h5ads) > 0
  Int n_samples = length(gex_h5ads)
  Boolean inputs_paired = n_samples > 0 && length(input_ids) == n_samples && (!multiome || length(atac_h5ads) == n_samples)
  Boolean can_label = defined(scanvi_model) || defined(ref_h5ad)

  if (!inputs_paired) {
    call utils.ErrorWithMessage as ErrorInputsNotPaired {
        input:
          message = "input_ids (~{length(input_ids)}) and gex_h5ads (~{n_samples}) must be non-empty and the same length, and atac_h5ads (~{length(atac_h5ads)}) must be empty (GEX-only) or the same length (multiome). On Terra, set members with no atac_h5ad value shorten atac_h5ads, so run multiome and GEX-only samples as separate sets."
    }
  }
  if (!can_label) {
    call utils.ErrorWithMessage as ErrorNoModelOrReference {
        input:
          message = "Give scanvi_model (to label with a saved model) or ref_h5ad (to train on the reference)."
    }
  }

  # Disk for the label-transfer machine: the inputs, each sample's outputs and scratch, and, when
  # training, the reference plus one reference-sized predictions file per sample.
  Float sample_gib = size(gex_h5ads, "GiB") + size(atac_h5ads, "GiB")
  Float ref_gib = if defined(scanvi_model) then 0.0 else size(ref_h5ad, "GiB")
  Int label_transfer_disk = select_first([disk_size, ceil(4 * sample_gib + (3 + n_samples) * ref_gib + size(scanvi_model, "GiB")) + 50])

  if (inputs_paired && can_label) {
    # Step 1 (training only): prepare the reference once for every sample
    if (!defined(scanvi_model)) {
      call PreprocessReference {
          input:
            ref_h5ad = select_first([ref_h5ad]),
            ref_label_column = ref_label_column,
            ref_batch_column = ref_batch_column,
            docker = docker
      }
    }

    # Step 2: preprocess and label every sample, one after another, on one machine. The GPU and
    # CPU tasks run the SAME command; they differ ONLY in the runtime GPU attributes, which
    # Cromwell cannot conditionally omit from a single task, so route by gpu_count.
    if (gpu_count > 0) {
      call MultiomeLabelTransfer {
          input:
            input_ids = input_ids,
            gex_h5ads = gex_h5ads,
            atac_h5ads = atac_h5ads,
            genome = genome,
            scanvi_model = scanvi_model,
            ref_h5ad = PreprocessReference.preprocessed_ref_h5ad,
            max_epochs = max_epochs,
            query_max_epochs = query_max_epochs,
            batch_size = batch_size,
            seed = seed,
            output_max_probability = output_max_probability,
            gpu_count = gpu_count,
            mem_size = mem_size,
            nthreads = nthreads,
            disk_size = label_transfer_disk,
            docker = docker
      }
    }
    if (gpu_count == 0) {
      call MultiomeLabelTransferCpu {
          input:
            input_ids = input_ids,
            gex_h5ads = gex_h5ads,
            atac_h5ads = atac_h5ads,
            genome = genome,
            scanvi_model = scanvi_model,
            ref_h5ad = PreprocessReference.preprocessed_ref_h5ad,
            max_epochs = max_epochs,
            query_max_epochs = query_max_epochs,
            batch_size = batch_size,
            seed = seed,
            output_max_probability = output_max_probability,
            mem_size = mem_size,
            nthreads = nthreads,
            disk_size = label_transfer_disk,
            docker = docker
      }
    }
  }

  output {
      # One entry per sample, in input order. atac_annotated_h5ads is empty in GEX-only mode.
      # Exactly one of the two label-transfer tasks runs; select_first picks it.
      Array[File] scanvi_predictions_h5ads = select_first([MultiomeLabelTransfer.scanvi_predictions_h5ads, MultiomeLabelTransferCpu.scanvi_predictions_h5ads])
      Array[File] gex_annotated_h5ads = select_first([MultiomeLabelTransfer.gex_annotated_h5ads, MultiomeLabelTransferCpu.gex_annotated_h5ads])
      Array[File] atac_annotated_h5ads = select_first([MultiomeLabelTransfer.atac_annotated_h5ads, MultiomeLabelTransferCpu.atac_annotated_h5ads])
      Array[File] scanvi_models_out = select_first([MultiomeLabelTransfer.scanvi_models_out, MultiomeLabelTransferCpu.scanvi_models_out])
      String pipeline_version_out = pipeline_version
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task 1: PreprocessReference (CPU-only, once per run, training only)
#
# Adapts the annotated reference into the form the trainer expects:
#   - Materializes counts from .raw for AIT-schema references (which have no .X)
#   - Maps the chosen label/batch obs columns to final_annotation/batch
#   - Tags the reference cells with modality rna_annotated
#   - Records the reference's id (its title or file name, plus a sha256 prefix of the file) in
#     uns['scanvi_reference']; every model trained on it carries that id as its name
# A run with a supplied scanvi_model skips it: the model needs no reference.
# ──────────────────────────────────────────────────────────────────────────────
task PreprocessReference {
    input {
        File ref_h5ad

        # Reference adaptation (AIT-schema support)
        String? ref_label_column
        String? ref_batch_column

        # Runtime attributes
        String docker
        Int disk_size = ceil(3 * size(ref_h5ad, "GiB")) + 20
        Int mem_size = 104
        Int nthreads = 16
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
import sys
import anndata as ad
import scanpy as sc

sys.path.insert(0, "/usr/local")
from reference_id import reference_id

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
ref.uns["scanvi_reference"] = reference_id("~{ref_h5ad}", ref.uns)
print(f"  Reference id: {ref.uns['scanvi_reference']}")

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
# Task 2: MultiomeLabelTransfer — every sample, one after another, on one machine.
#
# Checks every input first, then runs scanvi_sample.py once per sample: preprocess the sample's
# GEX (and ATAC), obtain the SCANVI model (supplied, mapped, or trained), transfer labels, and
# write the SCP-ready annotated matrices, the predictions and the model. Each output records the
# sample's id, input file name and model file name in uns. No sample metadata is added.
# Two variants that are IDENTICAL except for the runtime GPU attributes (Cromwell cannot
# conditionally omit them from a single task):
#   - MultiomeLabelTransfer     requests GPUs (gpu_count > 0)
#   - MultiomeLabelTransferCpu  is CPU-only (gpu_count = 0)
# scvi-tools auto-detects the accelerator. KEEP the two tasks' input / command / output blocks IN SYNC.
# ──────────────────────────────────────────────────────────────────────────────
task MultiomeLabelTransfer {
    input {
        # Every sample of the run, paired by position; atac_h5ads is empty in GEX-only mode
        Array[String] input_ids
        Array[File] gex_h5ads
        Array[File] atac_h5ads = []
        String genome = "hg38"

        # A supplied model, or the preprocessed reference to train on
        File? scanvi_model
        File? ref_h5ad

        Int? max_epochs
        Int? query_max_epochs
        Int batch_size = 128
        Int? seed
        Boolean output_max_probability = false

        # Runtime attributes
        String docker
        Int gpu_count = 2
        Int disk_size = 100
        Int mem_size = 104
        Int nthreads = 16
    }

    parameter_meta {
        input_ids: "One id per sample; it prefixes the sample's output file names and is recorded in each output's uns."
        gex_h5ads: "The samples' GEX h5ads (Optimus output), in the same order as input_ids."
        atac_h5ads: "The samples' ATAC cell-by-bin h5ads for multiome mode, one per GEX h5ad; empty for GEX-only mode."
        genome: "Genome for the ATAC gene-activity conversion: hg38, mm10 or mm39."
        scanvi_model: "Saved SCANVI model (.tar.gz) used for every sample instead of training."
        ref_h5ad: "Preprocessed reference, needed only when no scanvi_model is given."
        max_epochs: "SCVI/SCANVI training epochs when training (container default 500)."
        query_max_epochs: "Epochs when mapping a sample onto a reference-only model (container default 100)."
        batch_size: "SCVI/SCANVI minibatch size."
        seed: "Random seed for reproducible training; unset keeps training nondeterministic."
        output_max_probability: "Also write a max_probability obs column (the assigned label's confidence)."
        docker: "Docker image containing the scvi-scanvi runtime environment."
        gpu_count: "GPUs to request (default 2)."
        disk_size: "Disk size in GB."
        mem_size: "Memory in GiB; it only has to fit the largest sample, because samples run one at a time."
        nthreads: "CPUs."
    }

    command <<<
        set -euo pipefail
        mapfile -t IDS < ~{write_lines(input_ids)}
        mapfile -t GEX < ~{write_lines(gex_h5ads)}
        mapfile -t ATAC < ~{write_lines(atac_h5ads)}

        # Fail in seconds, naming the sample, if any input is missing, unreadable or mismatched.
        python3 /usr/local/validate_inputs.py --input-ids "${IDS[@]}" --gex "${GEX[@]}" --atac "${ATAC[@]}"

        # One sample at a time, each in its own process and scratch directory, so memory and disk
        # are freed before the next. Outputs are listed in input order, so outputs[i] is input i's.
        OUT=$(pwd)
        : > predictions.txt; : > gex_annotated.txt; : > atac_annotated.txt; : > models.txt
        for i in "${!IDS[@]}"; do
            id="${IDS[$i]}"
            atac_args=()
            if [ "${#ATAC[@]}" -gt 0 ]; then atac_args=(--atac "${ATAC[$i]}"); fi

            mkdir "work_$i" && cd "work_$i"
            ln -s /usr/local/gencode.v41.basic.annotation.gff3.gz .
            python3 /usr/local/scanvi_sample.py \
                --input-id "$id" \
                --gex "${GEX[$i]}" \
                "${atac_args[@]}" \
                --genome ~{genome} \
                ~{"--scanvi-model " + scanvi_model} \
                ~{"--ref " + ref_h5ad} \
                ~{"--max-epochs " + max_epochs} \
                ~{"--query-max-epochs " + query_max_epochs} \
                --batch-size ~{batch_size} \
                ~{"--seed " + seed} \
                ~{true="--output-max-probability" false="" output_max_probability}
            mv "${id}"_* "$OUT"/
            cd "$OUT" && rm -rf "work_$i"

            echo "${id}_SCANVI_predictions.h5ad" >> predictions.txt
            echo "${id}_gex_annotated_matrix.h5ad" >> gex_annotated.txt
            echo "${id}_scanvi_model.tar.gz" >> models.txt
            if [ "${#ATAC[@]}" -gt 0 ]; then echo "${id}_atac_annotated_matrix.h5ad" >> atac_annotated.txt; fi
        done
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
        # One retry. On Terra, "Retry with more memory" turns a retry after an out-of-memory kill
        # into a run with more memory; no other change is needed here.
        maxRetries: 1
    }

    output {
        Array[File] scanvi_predictions_h5ads = read_lines("predictions.txt")
        Array[File] gex_annotated_h5ads      = read_lines("gex_annotated.txt")
        Array[File] atac_annotated_h5ads     = read_lines("atac_annotated.txt")
        Array[File] scanvi_models_out        = read_lines("models.txt")
    }
}


# CPU-only variant of MultiomeLabelTransfer: identical inputs/command/outputs, but NO GPU runtime
# attributes and no gpu_count input. Keep in sync with MultiomeLabelTransfer above.
task MultiomeLabelTransferCpu {
    input {
        # Every sample of the run, paired by position; atac_h5ads is empty in GEX-only mode
        Array[String] input_ids
        Array[File] gex_h5ads
        Array[File] atac_h5ads = []
        String genome = "hg38"

        # A supplied model, or the preprocessed reference to train on
        File? scanvi_model
        File? ref_h5ad

        Int? max_epochs
        Int? query_max_epochs
        Int batch_size = 128
        Int? seed
        Boolean output_max_probability = false

        # Runtime attributes
        String docker
        Int disk_size = 100
        Int mem_size = 104
        Int nthreads = 16
    }

    parameter_meta {
        input_ids: "One id per sample; it prefixes the sample's output file names and is recorded in each output's uns."
        gex_h5ads: "The samples' GEX h5ads (Optimus output), in the same order as input_ids."
        atac_h5ads: "The samples' ATAC cell-by-bin h5ads for multiome mode, one per GEX h5ad; empty for GEX-only mode."
        genome: "Genome for the ATAC gene-activity conversion: hg38, mm10 or mm39."
        scanvi_model: "Saved SCANVI model (.tar.gz) used for every sample instead of training."
        ref_h5ad: "Preprocessed reference, needed only when no scanvi_model is given."
        max_epochs: "SCVI/SCANVI training epochs when training (container default 500)."
        query_max_epochs: "Epochs when mapping a sample onto a reference-only model (container default 100)."
        batch_size: "SCVI/SCANVI minibatch size."
        seed: "Random seed for reproducible training; unset keeps training nondeterministic."
        output_max_probability: "Also write a max_probability obs column (the assigned label's confidence)."
        docker: "Docker image containing the scvi-scanvi runtime environment."
        disk_size: "Disk size in GB."
        mem_size: "Memory in GiB; it only has to fit the largest sample, because samples run one at a time."
        nthreads: "CPUs."
    }

    command <<<
        set -euo pipefail
        mapfile -t IDS < ~{write_lines(input_ids)}
        mapfile -t GEX < ~{write_lines(gex_h5ads)}
        mapfile -t ATAC < ~{write_lines(atac_h5ads)}

        # Fail in seconds, naming the sample, if any input is missing, unreadable or mismatched.
        python3 /usr/local/validate_inputs.py --input-ids "${IDS[@]}" --gex "${GEX[@]}" --atac "${ATAC[@]}"

        # One sample at a time, each in its own process and scratch directory, so memory and disk
        # are freed before the next. Outputs are listed in input order, so outputs[i] is input i's.
        OUT=$(pwd)
        : > predictions.txt; : > gex_annotated.txt; : > atac_annotated.txt; : > models.txt
        for i in "${!IDS[@]}"; do
            id="${IDS[$i]}"
            atac_args=()
            if [ "${#ATAC[@]}" -gt 0 ]; then atac_args=(--atac "${ATAC[$i]}"); fi

            mkdir "work_$i" && cd "work_$i"
            ln -s /usr/local/gencode.v41.basic.annotation.gff3.gz .
            python3 /usr/local/scanvi_sample.py \
                --input-id "$id" \
                --gex "${GEX[$i]}" \
                "${atac_args[@]}" \
                --genome ~{genome} \
                ~{"--scanvi-model " + scanvi_model} \
                ~{"--ref " + ref_h5ad} \
                ~{"--max-epochs " + max_epochs} \
                ~{"--query-max-epochs " + query_max_epochs} \
                --batch-size ~{batch_size} \
                ~{"--seed " + seed} \
                ~{true="--output-max-probability" false="" output_max_probability}
            mv "${id}"_* "$OUT"/
            cd "$OUT" && rm -rf "work_$i"

            echo "${id}_SCANVI_predictions.h5ad" >> predictions.txt
            echo "${id}_gex_annotated_matrix.h5ad" >> gex_annotated.txt
            echo "${id}_scanvi_model.tar.gz" >> models.txt
            if [ "${#ATAC[@]}" -gt 0 ]; then echo "${id}_atac_annotated_matrix.h5ad" >> atac_annotated.txt; fi
        done
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        memory: "${mem_size} GiB"
        cpu: nthreads
        # One retry. On Terra, "Retry with more memory" turns a retry after an out-of-memory kill
        # into a run with more memory; no other change is needed here.
        maxRetries: 1
    }

    output {
        Array[File] scanvi_predictions_h5ads = read_lines("predictions.txt")
        Array[File] gex_annotated_h5ads      = read_lines("gex_annotated.txt")
        Array[File] atac_annotated_h5ads     = read_lines("atac_annotated.txt")
        Array[File] scanvi_models_out        = read_lines("models.txt")
    }
}
