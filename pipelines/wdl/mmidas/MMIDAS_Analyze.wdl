version 1.0

workflow MMIDAS_Analyze {

  meta {
    description: "Stage 3b-5 of the MMIDAS pipeline. Takes the confirmed evaluation_results.json from MMIDAS_Train (after human review of model_order), runs RF classification and silhouette analysis (03b) and state traversal preparation (03c) in parallel, then generates clusterability figures (04) and state traversal figures (05)."
    allowNestedInputs: true
  }

  String pipeline_version = "1.0.0"

  input {
    # ── Inputs from MMIDAS_Train (after human review) ─────────────────────────
    File   preprocessed_h5ad        # from MMIDAS_DataPrep
    File   checkpoints_manifest     # from MMIDAS_Train: TrainMixVAE output
    File   model_tar                # from MMIDAS_Train: TrainMixVAE output
    File   evaluation_results_json  # from MMIDAS_Train: Evaluate output — REVIEW THIS

    # ── Optional reference files ──────────────────────────────────────────────
    File?  kegg_toml    # KEGG/KEGG.toml — enables pathway box plots in 05
    File?  htree_file   # hierarchical taxonomy tree CSV — enables tree ordering in 03c

    # ── 03b: Classification parameters ───────────────────────────────────────
    Int    n_pca        = 100
    Int    k_fold       = 10
    String date_tag     = ""   # defaults to today's date inside the script

    # ── 03c: Traversal parameters ─────────────────────────────────────────────
    Int    n_traversal_steps  = 50
    Float  traversal_std_range = 3.0

    # ── 05: State traversal figure parameters ─────────────────────────────────
    Int    traversal_arm          = 0
    Int    n_selected_cats        = 10   # 0 = all active categories

    # ── Shared ────────────────────────────────────────────────────────────────
    Int    batch_size   = 5000
    Int    seed         = 0

    # ── Runtime ──────────────────────────────────────────────────────────────
    String docker           = "us.gcr.io/broad-gotc-prod/mmidas:1.0.0-0.1.0-1783366239"
    Int    analyze_disk_size = 200
    Int    analyze_mem_size  = 64
    Int    analyze_cpu       = 8
    Int    fig_disk_size     = 100
    Int    fig_mem_size      = 32
    Int    fig_cpu           = 4
  }

  # ── Restore model checkpoints (shared setup for 03b and 03c) ─────────────
  # Both 03b and 03c need the model/ directory. We untar once here so the
  # localized directory is available to both parallel tasks via their inputs.
  call RestoreModel {
    input:
      checkpoints_manifest    = checkpoints_manifest,
      model_tar               = model_tar,
      evaluation_results_json = evaluation_results_json,
      docker                  = docker,
      disk_size               = analyze_disk_size,
      mem_size                = 16,
      cpu                     = 2
  }

  # ── 03b and 03c run in parallel ────────────────────────────────────────────
  call Classify {
    input:
      anndata_h5ad            = preprocessed_h5ad,
      checkpoints_manifest    = RestoreModel.patched_manifest,
      model_dir_tar           = RestoreModel.model_dir_tar,
      evaluation_results_json = RestoreModel.patched_evaluation_results,
      n_pca                   = n_pca,
      k_fold                  = k_fold,
      batch_size              = batch_size,
      seed                    = seed,
      date_tag                = date_tag,
      docker                  = docker,
      disk_size               = analyze_disk_size,
      mem_size                = analyze_mem_size,
      cpu                     = analyze_cpu
  }

  call TraversalPrep {
    input:
      anndata_h5ad            = preprocessed_h5ad,
      checkpoints_manifest    = RestoreModel.patched_manifest,
      model_dir_tar           = RestoreModel.model_dir_tar,
      evaluation_results_json = RestoreModel.patched_evaluation_results,
      kegg_toml               = kegg_toml,
      htree_file              = htree_file,
      n_traversal_steps       = n_traversal_steps,
      traversal_std_range     = traversal_std_range,
      batch_size              = batch_size,
      seed                    = seed,
      docker                  = docker,
      disk_size               = analyze_disk_size,
      mem_size                = analyze_mem_size,
      cpu                     = analyze_cpu
  }

  # ── 04: Clusterability figures (depends on 03b) ────────────────────────────
  call Clusterability {
    input:
      classify_manifest = Classify.classify_manifest,
      clustering_tar    = Classify.clustering_tar,
      docker            = docker,
      disk_size         = fig_disk_size,
      mem_size          = fig_mem_size,
      cpu               = fig_cpu
  }

  # ── 05: State traversal figures (depends on 03c) ──────────────────────────
  call StateTraversal {
    input:
      anndata_h5ad            = preprocessed_h5ad,
      checkpoints_manifest    = RestoreModel.patched_manifest,
      model_dir_tar           = RestoreModel.model_dir_tar,
      evaluation_results_json = RestoreModel.patched_evaluation_results,
      traversal_manifest      = TraversalPrep.traversal_manifest,
      traversal_tar           = TraversalPrep.traversal_tar,
      arm                     = traversal_arm,
      n_selected_cats         = n_selected_cats,
      batch_size              = batch_size,
      seed                    = seed,
      docker                  = docker,
      disk_size               = fig_disk_size,
      mem_size                = fig_mem_size,
      cpu                     = fig_cpu
  }

  output {
    # Clusterability (04)
    Array[File] clusterability_figures     = Clusterability.figures
    File        clusterability_manifest    = Clusterability.clusterability_manifest

    # State traversal (05)
    Array[File] state_traversal_figures    = StateTraversal.figures
    File        state_traversal_manifest   = StateTraversal.state_traversal_manifest

    String      pipeline_version_out       = pipeline_version
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task: RestoreModel
#
# Unpacks model_tar and rewrites absolute paths in checkpoints_manifest.json
# to point to the current working directory. The patched manifest and a
# re-tarred model directory are passed to Classify, TraversalPrep, and
# StateTraversal so each task can reconstruct the model locally.
# ──────────────────────────────────────────────────────────────────────────────
task RestoreModel {
  input {
    File   checkpoints_manifest
    File   model_tar
    File   evaluation_results_json
    String docker
    Int    disk_size
    Int    mem_size
    Int    cpu
  }

  command <<<
    set -euo pipefail
    mkdir -p out/model

    tar -xzf "~{model_tar}" -C out/

    ABS_OUT=$(realpath out)

  # Normalize checkpoints_manifest.json to the local execution directory.
  # Be resilient to legacy/new schema differences (e.g. missing output_dir).
  python3 - <<CODE
import glob
import json
import os
import sys

src = "~{checkpoints_manifest}"
dst = "out/checkpoints_manifest.json"
abs_out = os.path.realpath("out")
model_dir = os.path.join(abs_out, "model")

with open(src) as fh:
  manifest = json.load(fh)

local_ckpts = sorted(glob.glob(os.path.join(model_dir, "*.pth")))
if not local_ckpts:
  print("ERROR: no .pth checkpoints found under out/model after extracting model_tar", file=sys.stderr)
  sys.exit(1)

manifest["output_dir"] = abs_out
manifest["checkpoints"] = local_ckpts

with open(dst, "w") as fh:
  json.dump(manifest, fh, indent=2)
CODE

  # Normalize evaluation_results.json selected_model/output_dir paths.
  python3 - <<CODE
import glob
import json
import os

src = "~{evaluation_results_json}"
dst = "out/evaluation_results.json"
abs_out = os.path.realpath("out")
model_dir = os.path.join(abs_out, "model")
local_ckpts = sorted(glob.glob(os.path.join(model_dir, "*.pth")))

with open(src) as fh:
  eval_results = json.load(fh)

eval_results["output_dir"] = abs_out

selected = eval_results.get("selected_model")
if selected:
  selected_base = os.path.basename(selected)
  selected_local = os.path.join(model_dir, selected_base)
  if os.path.exists(selected_local):
    eval_results["selected_model"] = selected_local
  elif local_ckpts:
    # Fallback to last checkpoint if exact basename is unavailable.
    eval_results["selected_model"] = local_ckpts[-1]

with open(dst, "w") as fh:
  json.dump(eval_results, fh, indent=2)
CODE

    # Re-tar so downstream tasks receive a self-contained model bundle
    tar -czf model_dir.tar.gz -C out model/
  >>>

  output {
    File patched_manifest           = "out/checkpoints_manifest.json"
    File patched_evaluation_results = "out/evaluation_results.json"
    File model_dir_tar              = "model_dir.tar.gz"
  }

  runtime {
    docker:      docker
    disks:       "local-disk ~{disk_size} HDD"
    memory:      "~{mem_size} GiB"
    cpu:         cpu
    preemptible: 3
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task: Classify  (03b)
#
# RF classification + silhouette analysis across label sets and embeddings.
# Outputs clustering pickles (tarred) and classify_manifest.json.
# ──────────────────────────────────────────────────────────────────────────────
task Classify {
  input {
    File   anndata_h5ad
    File   checkpoints_manifest
    File   model_dir_tar
    File   evaluation_results_json
    Int    n_pca
    Int    k_fold
    Int    batch_size
    Int    seed
    String date_tag
    String docker
    Int    disk_size
    Int    mem_size
    Int    cpu
  }

  parameter_meta {
    anndata_h5ad:            "Preprocessed AnnData .h5ad file."
    checkpoints_manifest:    "Patched checkpoints_manifest.json from RestoreModel."
    model_dir_tar:           "Tarball of the model/ checkpoint directory from RestoreModel."
    evaluation_results_json: "evaluation_results.json from MMIDAS_Train Evaluate task."
    n_pca:                   "Number of PCA components for the linear embedding baseline."
    k_fold:                  "Number of folds for k-fold cross-validation."
    batch_size:              "Batch size for model inference."
    seed:                    "Random seed (must match training)."
    date_tag:                "Date string embedded in pickle filenames (empty = today)."
    docker:                  "MMIDAS Docker image URI."
    disk_size:               "Boot disk size in GB."
    mem_size:                "Memory in GB."
    cpu:                     "Number of CPU cores."
  }

  command <<<
    set -euo pipefail
    mkdir -p out

    tar -xzf "~{model_dir_tar}" -C out/

    python3 /usr/local/03b_classify.py \
      --anndata_h5ad            "~{anndata_h5ad}" \
      --checkpoints_manifest    "~{checkpoints_manifest}" \
      --evaluation_results_json "~{evaluation_results_json}" \
      --output_dir              out \
      --n_pca                   ~{n_pca} \
      --k_fold                  ~{k_fold} \
      --batch_size              ~{batch_size} \
      --seed                    ~{seed} \
      ~{if date_tag != "" then "--date_tag \"" + date_tag + "\"" else ""}

    tar -czf clustering.tar.gz -C out clustering/
  >>>

  output {
    File classify_manifest = "out/classify_manifest.json"
    File clustering_tar    = "clustering.tar.gz"
  }

  runtime {
    docker:      docker
    disks:       "local-disk ~{disk_size} HDD"
    memory:      "~{mem_size} GiB"
    cpu:         cpu
    preemptible: 1
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task: TraversalPrep  (03c)
#
# Pre-computes state traversal data: per-gene variation scores, KEGG pathway
# gene-index mapping, taxonomy ordering, and state-space traversal paths.
# Outputs traversal pickles (tarred) and traversal_manifest.json.
# ──────────────────────────────────────────────────────────────────────────────
task TraversalPrep {
  input {
    File   anndata_h5ad
    File   checkpoints_manifest
    File   model_dir_tar
    File   evaluation_results_json
    File?  kegg_toml
    File?  htree_file
    Int    n_traversal_steps
    Float  traversal_std_range
    Int    batch_size
    Int    seed
    String docker
    Int    disk_size
    Int    mem_size
    Int    cpu
  }

  parameter_meta {
    anndata_h5ad:            "Preprocessed AnnData .h5ad file."
    checkpoints_manifest:    "Patched checkpoints_manifest.json from RestoreModel."
    model_dir_tar:           "Tarball of the model/ checkpoint directory from RestoreModel."
    evaluation_results_json: "evaluation_results.json from MMIDAS_Train Evaluate task."
    kegg_toml:               "Optional KEGG.toml for pathway gene-set mapping."
    htree_file:              "Optional hierarchical taxonomy tree CSV for category ordering."
    n_traversal_steps:       "Number of points along the state traversal path."
    traversal_std_range:     "Traversal spans ± (this value) × std along PC1."
    batch_size:              "Batch size for model inference."
    seed:                    "Random seed."
    docker:                  "MMIDAS Docker image URI."
    disk_size:               "Boot disk size in GB."
    mem_size:                "Memory in GB."
    cpu:                     "Number of CPU cores."
  }

  command <<<
    set -euo pipefail
    mkdir -p out/state

    tar -xzf "~{model_dir_tar}" -C out/

    python3 /usr/local/03c_traversal_prep.py \
      --anndata_h5ad            "~{anndata_h5ad}" \
      --checkpoints_manifest    "~{checkpoints_manifest}" \
      --evaluation_results_json "~{evaluation_results_json}" \
      --output_dir              out \
      ~{if defined(kegg_toml)  then "--kegg_toml \""  + select_first([kegg_toml])  + "\"" else ""} \
      ~{if defined(htree_file) then "--htree_file \"" + select_first([htree_file]) + "\"" else ""} \
      --n_traversal_steps       ~{n_traversal_steps} \
      --traversal_std_range     ~{traversal_std_range} \
      --batch_size              ~{batch_size} \
      --seed                    ~{seed}

    tar -czf traversal.tar.gz -C out state/
    cp out/taxonomy_order_K_*.npy .
  >>>

  output {
    File        traversal_manifest = "out/traversal_manifest.json"
    File        traversal_tar      = "traversal.tar.gz"
    Array[File] taxonomy_order     = glob("taxonomy_order_K_*.npy")
  }

  runtime {
    docker:      docker
    disks:       "local-disk ~{disk_size} HDD"
    memory:      "~{mem_size} GiB"
    cpu:         cpu
    preemptible: 1
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task: Clusterability  (04)
#
# Generates classification accuracy bar chart, silhouette score curves, and
# confusion matrix heatmaps from the 03b clustering pickles.
# ──────────────────────────────────────────────────────────────────────────────
task Clusterability {
  input {
    File   classify_manifest
    File   clustering_tar
    String docker
    Int    disk_size
    Int    mem_size
    Int    cpu
  }

  parameter_meta {
    classify_manifest: "classify_manifest.json from Classify task."
    clustering_tar:    "Tarball of the clustering/ pickles directory from Classify task."
    docker:            "MMIDAS Docker image URI."
    disk_size:         "Boot disk size in GB."
    mem_size:          "Memory in GB."
    cpu:               "Number of CPU cores."
  }

  command <<<
    set -euo pipefail
    mkdir -p out

    # Restore clustering pickles and rewrite the manifest path
    tar -xzf "~{clustering_tar}" -C out/

    cp "~{classify_manifest}" out/classify_manifest.json
    ORIG_CLUST=$(python3 -c "import json; print(json.load(open('out/classify_manifest.json'))['clustering_dir'])")
    ABS_OUT=$(realpath out)
    sed -i "s|${ORIG_CLUST}|${ABS_OUT}/clustering|g" out/classify_manifest.json
    sed -i "s|\"output_dir\": \"[^\"]*\"|\"output_dir\": \"${ABS_OUT}\"|" out/classify_manifest.json

    python3 /usr/local/04_clusterability.py \
      --classify_manifest "out/classify_manifest.json" \
      --output_dir        out
  >>>

  output {
    Array[File] figures                = glob("out/*.png")
    File        clusterability_manifest = "out/clusterability_manifest.json"
  }

  runtime {
    docker:      docker
    disks:       "local-disk ~{disk_size} HDD"
    memory:      "~{mem_size} GiB"
    cpu:         cpu
    preemptible: 3
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task: StateTraversal  (05)
#
# Generates per-category state-space scatter plots and, when KEGG pathways
# were provided to TraversalPrep, per-pathway and per-category box plots.
# ──────────────────────────────────────────────────────────────────────────────
task StateTraversal {
  input {
    File   anndata_h5ad
    File   checkpoints_manifest
    File   model_dir_tar
    File   evaluation_results_json
    File   traversal_manifest
    File   traversal_tar
    Int    arm
    Int    n_selected_cats
    Int    batch_size
    Int    seed
    String docker
    Int    disk_size
    Int    mem_size
    Int    cpu
  }

  parameter_meta {
    anndata_h5ad:            "Preprocessed AnnData .h5ad file."
    checkpoints_manifest:    "Patched checkpoints_manifest.json from RestoreModel."
    model_dir_tar:           "Tarball of the model/ checkpoint directory from RestoreModel."
    evaluation_results_json: "evaluation_results.json from MMIDAS_Train Evaluate task."
    traversal_manifest:      "traversal_manifest.json from TraversalPrep task."
    traversal_tar:           "Tarball of the state/ traversal pickle directory from TraversalPrep."
    arm:                     "Which encoder arm to use for state-space scatter plots."
    n_selected_cats:         "Number of active categories to generate figures for (0 = all)."
    batch_size:              "Batch size for model inference."
    seed:                    "Random seed."
    docker:                  "MMIDAS Docker image URI."
    disk_size:               "Boot disk size in GB."
    mem_size:                "Memory in GB."
    cpu:                     "Number of CPU cores."
  }

  command <<<
    set -euo pipefail
    mkdir -p out/state

    tar -xzf "~{model_dir_tar}"  -C out/
    tar -xzf "~{traversal_tar}"  -C out/

    # Rewrite absolute paths in both manifests
    cp "~{checkpoints_manifest}" out/checkpoints_manifest.json
    cp "~{traversal_manifest}"   out/traversal_manifest.json
    ABS_OUT=$(realpath out)

    ORIG_CKPT=$(python3 -c "import json; print(json.load(open('out/checkpoints_manifest.json'))['output_dir'])")
    sed -i "s|${ORIG_CKPT}|${ABS_OUT}|g" out/checkpoints_manifest.json

    ORIG_TRAV=$(python3 -c "import json; print(json.load(open('out/traversal_manifest.json'))['output_dir'])")
    sed -i "s|${ORIG_TRAV}|${ABS_OUT}|g" out/traversal_manifest.json

    python3 /usr/local/05_state_traversal.py \
      --anndata_h5ad            "~{anndata_h5ad}" \
      --checkpoints_manifest    "out/checkpoints_manifest.json" \
      --evaluation_results_json "~{evaluation_results_json}" \
      --traversal_manifest      "out/traversal_manifest.json" \
      --arm                     ~{arm} \
      --n_selected_cats         ~{n_selected_cats} \
      --batch_size              ~{batch_size} \
      --seed                    ~{seed}

    # Collect all figures into the task root for glob output
    find out/state -name "*.png" -exec cp {} . \;
  >>>

  output {
    Array[File] figures                 = glob("*.png")
    File        state_traversal_manifest = "out/state_traversal_manifest.json"
  }

  runtime {
    docker:      docker
    disks:       "local-disk ~{disk_size} HDD"
    memory:      "~{mem_size} GiB"
    cpu:         cpu
    preemptible: 3
  }
}
