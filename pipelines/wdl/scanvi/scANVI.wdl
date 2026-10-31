version 1.0

import "../../../tasks/wdl/Utilities.wdl" as utils

workflow scANVI {

  meta {
    description: "Pipeline for cell type label transfer using SCVI and SCANVI models. Labels one or more samples, each independently and with its own outputs, at one or more label levels (e.g. class and subclass), either with a supplied model archive or by training on an annotated reference plus the sample. Each sample's single-cell RNA (GEX) and, optionally, ATAC data are labelled; each input file gives an <input file name>_annotated.h5ad written ready for Single Cell Portal. All samples of a run go through one machine, one after another."
    allowNestedInputs: true
  }

  input {
      # Samples, paired by position: input_ids[i] is the id of gex_h5ads[i] (and atac_h5ads[i]),
      # e.g. the library's id, and is recorded in each of its outputs' uns. Outputs are named
      # after the input files (<name>.h5ad -> <name>_annotated.h5ad), so the input files must
      # have distinct names. A single sample is a one-element array.
      Array[String] input_ids
      Array[File] gex_h5ads

      # ATAC cell-by-bin h5ads for multiome mode: one per gex_h5ads entry, in the same order.
      # Empty (default) runs every sample in GEX-only mode.
      Array[File] atac_h5ads = []

      # Optional model archive (.tar.gz: a reference's scVI model and one SCANVI model per label
      # level, saved without data), used for every sample; every level labels every sample. A
      # model this pipeline trained loads as is; a model trained on the reference alone
      # (warp-scanvi-training) gets each sample mapped onto it (scArches). Without a model, each
      # sample trains its own models on ref_h5ad plus that sample.
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

      # Reference adaptation when training. ref_label_columns are the label levels: the obs
      # columns to label each sample with, one SCANVI model each, e.g. ["class", "subclass"]; the
      # first gives the default UMAP. ref_batch_column is the batch. When unset, AIT references
      # default to subclass/donor_id and PBMC-style references to final_annotation/batch.
      Array[String] ref_label_columns = []
      String? ref_batch_column

      # Genome for the ATAC cell-by-bin -> gene-activity conversion (multiome mode only):
      # "hg38" (default, human), "mm10", or "mm39" (mouse).
      String genome = "hg38"

      # Each level's per-cell maximum SCANVI posterior probability (the confidence of the
      # assigned label), as a <level>_max_probability obs column. Set false to leave it out.
      Boolean output_max_probability = true

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
  String docker = "us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:59206dd15b3fb8170b63c0da429b35ba597847ee30e4cf520b513ed9bebad22e"

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
  # training, the reference and its scratch.
  Float sample_gib = size(gex_h5ads, "GiB") + size(atac_h5ads, "GiB")
  Float ref_gib = if defined(scanvi_model) then 0.0 else size(ref_h5ad, "GiB")
  Int label_transfer_disk = select_first([disk_size, ceil(4 * sample_gib + 3 * ref_gib + size(scanvi_model, "GiB")) + 50])

  # Each sample's output file names, in input order, named after its input files. They are the
  # label-transfer task's outputs, so outputs[i] belongs to input_ids[i], and Cromwell knows every
  # output before the task runs (on Google Cloud it saves only outputs it can name in advance).
  # These loops only build strings: they call no task and start no machines. All samples are
  # labelled on the one machine below. A run that trains also saves each sample's models.
  scatter (gex in gex_h5ads) {
    String gex_name   = basename(gex, ".h5ad") + "_annotated.h5ad"
    String model_name = basename(gex, ".h5ad") + "_scanvi_model.tar.gz"
  }
  scatter (atac in atac_h5ads) {
    String atac_name = basename(atac, ".h5ad") + "_annotated.h5ad"
  }
  Array[String] model_names = if defined(scanvi_model) then [] else model_name

  if (inputs_paired && can_label) {
    # Step 1 (training only): prepare the reference once for every sample
    if (!defined(scanvi_model)) {
      call PreprocessReference {
          input:
            ref_h5ad = select_first([ref_h5ad]),
            ref_label_columns = ref_label_columns,
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
            gex_names = gex_name,
            atac_names = atac_name,
            model_names = model_names,
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
            gex_names = gex_name,
            atac_names = atac_name,
            model_names = model_names,
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
      # One entry per sample, in input order. atac_annotated_h5ads is empty in GEX-only mode, and
      # scanvi_models_out is empty when a model was supplied (only a run that trains saves models).
      # Exactly one of the two label-transfer tasks runs; select_first picks it.
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
#   - Keeps the chosen label columns (the levels, listed in uns['scanvi_levels']) as strings and
#     maps the batch column to batch
#   - Records the reference's id (its title or file name, plus a sha256 prefix of the file) in
#     uns['scanvi_reference']; every model trained on it carries that id as its name
# A run with a supplied scanvi_model skips it: the model needs no reference.
# ──────────────────────────────────────────────────────────────────────────────
task PreprocessReference {
    input {
        File ref_h5ad

        # Reference adaptation (AIT-schema support)
        Array[String] ref_label_columns = []
        String? ref_batch_column

        # Runtime attributes
        String docker
        Int disk_size = ceil(3 * size(ref_h5ad, "GiB")) + 20
        Int mem_size = 104
        Int nthreads = 16
    }

    parameter_meta {
        ref_h5ad: "Annotated reference AnnData h5ad file. PBMC-style references carry cell type labels in obs['final_annotation'] and a batch in obs['batch']; AIT-schema references (uns['schema_version']+uns['hierarchy']) are auto-detected and adapted (see ref_label_columns/ref_batch_column)."
        ref_label_columns: "Reference obs columns to label with, one SCANVI model each (the label levels). When empty, defaults to ['subclass'] for AIT references and ['final_annotation'] otherwise."
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
from model_archive import check_levels, reference_id

# Reference-adaptation controls (from WDL inputs; empty = use defaults).
ref_label_columns = "~{sep=' ' ref_label_columns}".split()
ref_batch_column = "~{default='' ref_batch_column}".strip()


def normalize_reference(ref, label_cols_override, batch_col_override):
    """Adapt the annotated reference into the form the trainer expects: counts in
    .X, each label level's column as strings, its levels in uns['scanvi_levels'], and
    a batch in obs['batch'].

    Supports AIT-schema atlases (Allen Institute Taxonomy: no .X — counts live in
    .raw; labels in a class/subclass/cluster_id hierarchy; batch in donor_id) and
    PBMC-style references (already have .X + final_annotation + batch).
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
    levels = label_cols_override or (["subclass"] if is_ait else ["final_annotation"])
    batch_col = batch_col_override or ("donor_id" if is_ait else "batch")
    check_levels(levels)
    for col, role in [(c, "label") for c in levels] + [(batch_col, "batch")]:
        if col not in ref.obs.columns:
            raise SystemExit(
                f"ERROR: reference {role} column '{col}' not found in obs "
                f"(available: {sorted(ref.obs.columns)}).")

    for level in levels:
        ref.obs[level] = ref.obs[level].astype(str)
    ref.obs["batch"] = ref.obs[batch_col].astype(str)
    ref.uns["scanvi_levels"] = levels
    print(f"  Reference ready: levels {levels}, batch='{batch_col}' -> batch, shape {ref.shape}")
    return ref


print("Loading reference...")
ref = sc.read_h5ad("~{ref_h5ad}")
print(f"  Reference loaded: {ref.shape}")
# Adapt AIT-schema (or other) references into the trainer's expected form.
ref = normalize_reference(ref, ref_label_columns, ref_batch_column)
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
# GEX (and ATAC), obtain each label level's SCANVI model (supplied, mapped, or trained), transfer
# labels, and write the SCP-ready annotated matrices (and, when training, the models). Each output
# records the sample's id, input file name, model name and levels in uns. No sample metadata.
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
        # Each sample's output file names, in input order (built by the workflow from the inputs);
        # model_names is empty when a model is supplied
        Array[String] gex_names
        Array[String] atac_names
        Array[String] model_names
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
        Boolean output_max_probability = true

        # Runtime attributes
        String docker
        Int gpu_count = 2
        Int disk_size = 100
        Int mem_size = 104
        Int nthreads = 16
    }

    parameter_meta {
        input_ids: "One id per sample, recorded in each of its outputs' uns."
        gex_names: "Output file names, one per sample in input order; set by the workflow."
        atac_names: "Output file names, one per sample in input order (empty in GEX-only mode); set by the workflow."
        model_names: "Output file names, one per sample in input order (empty when a model is supplied); set by the workflow."
        gex_h5ads: "The samples' GEX h5ads (Optimus output), in the same order as input_ids."
        atac_h5ads: "The samples' ATAC cell-by-bin h5ads for multiome mode, one per GEX h5ad; empty for GEX-only mode."
        genome: "Genome for the ATAC gene-activity conversion: hg38, mm10 or mm39."
        scanvi_model: "Model archive (.tar.gz) used for every sample instead of training."
        ref_h5ad: "Preprocessed reference, needed only when no scanvi_model is given."
        max_epochs: "SCVI/SCANVI training epochs when training (container default 500)."
        query_max_epochs: "Epochs when mapping a sample onto a reference-only model (container default 100)."
        batch_size: "SCVI/SCANVI minibatch size."
        seed: "Random seed for reproducible training; unset keeps training nondeterministic."
        output_max_probability: "Write a <level>_max_probability obs column (each assigned label's confidence); default true."
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
        mapfile -t GEX_OUT < ~{write_lines(gex_names)}
        # Cromwell writes an empty list as one blank line, so read these only when there are some.
        ATAC=() ATAC_OUT=() MODEL_OUT=()
        if [ ~{length(atac_h5ads)} -gt 0 ]; then
            mapfile -t ATAC < ~{write_lines(atac_h5ads)}
            mapfile -t ATAC_OUT < ~{write_lines(atac_names)}
        fi
        if [ ~{length(model_names)} -gt 0 ]; then mapfile -t MODEL_OUT < ~{write_lines(model_names)}; fi

        # Fail in seconds, naming the sample, if any input is missing, unreadable or mismatched, or
        # two inputs share a file name (their outputs would collide).
        python3 /usr/local/validate_inputs.py --input-ids "${IDS[@]}" --gex "${GEX[@]}" --atac "${ATAC[@]}" \
            --outputs "${GEX_OUT[@]}" "${ATAC_OUT[@]}"

        # One sample at a time, each in its own process and scratch directory, so memory and disk
        # are freed before the next. Outputs are written straight to the task directory.
        OUT=$(pwd)
        for i in "${!IDS[@]}"; do
            atac_args=() model_args=()
            if [ "${#ATAC[@]}" -gt 0 ]; then atac_args=(--atac "${ATAC[$i]}" --atac-out "$OUT/${ATAC_OUT[$i]}"); fi
            if [ "${#MODEL_OUT[@]}" -gt 0 ]; then model_args=(--model-out "$OUT/${MODEL_OUT[$i]}"); fi

            mkdir "work_$i" && cd "work_$i"
            python3 /usr/local/scanvi_sample.py \
                --input-id "${IDS[$i]}" \
                --gex "${GEX[$i]}" \
                --gex-out "$OUT/${GEX_OUT[$i]}" \
                "${atac_args[@]}" \
                "${model_args[@]}" \
                --genome ~{genome} \
                ~{"--scanvi-model " + scanvi_model} \
                ~{"--ref " + ref_h5ad} \
                ~{"--max-epochs " + max_epochs} \
                ~{"--query-max-epochs " + query_max_epochs} \
                --batch-size ~{batch_size} \
                ~{"--seed " + seed} \
                ~{true="--output-max-probability" false="" output_max_probability}
            cd "$OUT" && rm -rf "work_$i"
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
        # The names come from the workflow, in input order; a missing file fails the task.
        Array[File] gex_annotated_h5ads      = gex_names
        Array[File] atac_annotated_h5ads     = atac_names
        Array[File] scanvi_models_out        = model_names
    }
}


# CPU-only variant of MultiomeLabelTransfer: identical inputs/command/outputs, but NO GPU runtime
# attributes and no gpu_count input. Keep in sync with MultiomeLabelTransfer above.
task MultiomeLabelTransferCpu {
    input {
        # Every sample of the run, paired by position; atac_h5ads is empty in GEX-only mode
        Array[String] input_ids
        # Each sample's output file names, in input order (built by the workflow from the inputs);
        # model_names is empty when a model is supplied
        Array[String] gex_names
        Array[String] atac_names
        Array[String] model_names
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
        Boolean output_max_probability = true

        # Runtime attributes
        String docker
        Int disk_size = 100
        Int mem_size = 104
        Int nthreads = 16
    }

    parameter_meta {
        input_ids: "One id per sample, recorded in each of its outputs' uns."
        gex_names: "Output file names, one per sample in input order; set by the workflow."
        atac_names: "Output file names, one per sample in input order (empty in GEX-only mode); set by the workflow."
        model_names: "Output file names, one per sample in input order (empty when a model is supplied); set by the workflow."
        gex_h5ads: "The samples' GEX h5ads (Optimus output), in the same order as input_ids."
        atac_h5ads: "The samples' ATAC cell-by-bin h5ads for multiome mode, one per GEX h5ad; empty for GEX-only mode."
        genome: "Genome for the ATAC gene-activity conversion: hg38, mm10 or mm39."
        scanvi_model: "Model archive (.tar.gz) used for every sample instead of training."
        ref_h5ad: "Preprocessed reference, needed only when no scanvi_model is given."
        max_epochs: "SCVI/SCANVI training epochs when training (container default 500)."
        query_max_epochs: "Epochs when mapping a sample onto a reference-only model (container default 100)."
        batch_size: "SCVI/SCANVI minibatch size."
        seed: "Random seed for reproducible training; unset keeps training nondeterministic."
        output_max_probability: "Write a <level>_max_probability obs column (each assigned label's confidence); default true."
        docker: "Docker image containing the scvi-scanvi runtime environment."
        disk_size: "Disk size in GB."
        mem_size: "Memory in GiB; it only has to fit the largest sample, because samples run one at a time."
        nthreads: "CPUs."
    }

    command <<<
        set -euo pipefail
        mapfile -t IDS < ~{write_lines(input_ids)}
        mapfile -t GEX < ~{write_lines(gex_h5ads)}
        mapfile -t GEX_OUT < ~{write_lines(gex_names)}
        # Cromwell writes an empty list as one blank line, so read these only when there are some.
        ATAC=() ATAC_OUT=() MODEL_OUT=()
        if [ ~{length(atac_h5ads)} -gt 0 ]; then
            mapfile -t ATAC < ~{write_lines(atac_h5ads)}
            mapfile -t ATAC_OUT < ~{write_lines(atac_names)}
        fi
        if [ ~{length(model_names)} -gt 0 ]; then mapfile -t MODEL_OUT < ~{write_lines(model_names)}; fi

        # Fail in seconds, naming the sample, if any input is missing, unreadable or mismatched, or
        # two inputs share a file name (their outputs would collide).
        python3 /usr/local/validate_inputs.py --input-ids "${IDS[@]}" --gex "${GEX[@]}" --atac "${ATAC[@]}" \
            --outputs "${GEX_OUT[@]}" "${ATAC_OUT[@]}"

        # One sample at a time, each in its own process and scratch directory, so memory and disk
        # are freed before the next. Outputs are written straight to the task directory.
        OUT=$(pwd)
        for i in "${!IDS[@]}"; do
            atac_args=() model_args=()
            if [ "${#ATAC[@]}" -gt 0 ]; then atac_args=(--atac "${ATAC[$i]}" --atac-out "$OUT/${ATAC_OUT[$i]}"); fi
            if [ "${#MODEL_OUT[@]}" -gt 0 ]; then model_args=(--model-out "$OUT/${MODEL_OUT[$i]}"); fi

            mkdir "work_$i" && cd "work_$i"
            python3 /usr/local/scanvi_sample.py \
                --input-id "${IDS[$i]}" \
                --gex "${GEX[$i]}" \
                --gex-out "$OUT/${GEX_OUT[$i]}" \
                "${atac_args[@]}" \
                "${model_args[@]}" \
                --genome ~{genome} \
                ~{"--scanvi-model " + scanvi_model} \
                ~{"--ref " + ref_h5ad} \
                ~{"--max-epochs " + max_epochs} \
                ~{"--query-max-epochs " + query_max_epochs} \
                --batch-size ~{batch_size} \
                ~{"--seed " + seed} \
                ~{true="--output-max-probability" false="" output_max_probability}
            cd "$OUT" && rm -rf "work_$i"
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
        # The names come from the workflow, in input order; a missing file fails the task.
        Array[File] gex_annotated_h5ads      = gex_names
        Array[File] atac_annotated_h5ads     = atac_names
        Array[File] scanvi_models_out        = model_names
    }
}
