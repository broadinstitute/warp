version 1.0

workflow MMIDAS_Train {

  meta {
    description: "Stage 2+3a of the MMIDAS pipeline. Optionally trains a UDAGAN VAE-GAN augmenter (02a), trains the core cpl-mixVAE model with iterative pruning (02b), and evaluates all checkpoints to select the optimal model order (03a). Outputs evaluation figures and evaluation_results.json for human review before running MMIDAS_Analyze."
    allowNestedInputs: true
  }

  String pipeline_version = "1.0.0"

  input {
    # ── Input data ───────────────────────────────────────────────────────────
    File   preprocessed_h5ad      # output of MMIDAS_DataPrep

    # ── Augmenter (optional) ─────────────────────────────────────────────────
    Boolean run_augmenter          = false
    String  augmenter_tag          = "mmidas"
    Int     augmenter_z_dim        = 10
    Int     augmenter_noise_dim    = 50
    Int     augmenter_fc_dim       = 500
    Int     augmenter_n_epoch      = 500
    Int     augmenter_batch_size   = 512

    # ── cpl-mixVAE architecture ───────────────────────────────────────────────
    Int     n_categories           = 120
    Int     state_dim              = 2
    Int     n_arm                  = 2
    Int     latent_dim             = 10
    Int     fc_dim                 = 100
    Float   x_drop                 = 0.25
    Float   s_drop                 = 0.0
    Float   temp                   = 1.0
    Float   tau                    = 0.1
    Float   beta                   = 1.0
    Float   lam                    = 1.0
    Float   lam_pc                 = 1.0
    String  training_mode          = "MSE"   # MSE or ZINB

    # ── cpl-mixVAE training ───────────────────────────────────────────────────
    Int     n_epoch                = 10000
    Int     n_epoch_p              = 1000
    Float   min_con                = 0.99
    Int     max_prun_it            = 14
    Float   lr                     = 0.001
    Int     batch_size             = 5000
    Int     n_aug_smp              = 0
    Float   train_size             = 0.9
    Int     seed                   = 0

    # ── Evaluation ────────────────────────────────────────────────────────────
    Float   k_select_thr           = 0.95

    # ── Runtime ──────────────────────────────────────────────────────────────
    String  docker                 = "us.gcr.io/broad-gotc-prod/mmidas:1.0.0-0.1.0-1782763680"
    Int     train_disk_size        = 200
    Int     train_mem_size         = 64
    Int     train_cpu              = 8
    Int     eval_disk_size         = 100
    Int     eval_mem_size          = 32
    Int     eval_cpu               = 4
    # Set train_gpu to 1 to attach a GPU to the training task.
    # On Terra use: {"nvidia_tesla_t4": 1} in the runtime attributes.
    Int     train_gpu              = 0
  }

  # ── 02a: Optional augmenter training ──────────────────────────────────────
  if (run_augmenter) {
    call TrainAugmenter {
      input:
        anndata_h5ad = preprocessed_h5ad,
        tag          = augmenter_tag,
        z_dim        = augmenter_z_dim,
        noise_dim    = augmenter_noise_dim,
        fc_dim       = augmenter_fc_dim,
        n_epoch      = augmenter_n_epoch,
        batch_size   = augmenter_batch_size,
        use_gpu      = train_gpu > 0,
        docker       = docker,
        disk_size    = train_disk_size,
        mem_size     = train_mem_size,
        cpu          = train_cpu,
        gpu          = train_gpu
    }
  }

  # ── 02b: Core cpl-mixVAE training + pruning ────────────────────────────────
  call TrainMixVAE {
    input:
      anndata_h5ad          = preprocessed_h5ad,
      augmenter_checkpoint  = TrainAugmenter.augmenter_checkpoint,
      n_categories          = n_categories,
      state_dim             = state_dim,
      n_arm                 = n_arm,
      latent_dim            = latent_dim,
      fc_dim                = fc_dim,
      x_drop                = x_drop,
      s_drop                = s_drop,
      temp                  = temp,
      tau                   = tau,
      beta                  = beta,
      lam                   = lam,
      lam_pc                = lam_pc,
      training_mode         = training_mode,
      n_epoch               = n_epoch,
      n_epoch_p             = n_epoch_p,
      min_con               = min_con,
      max_prun_it           = max_prun_it,
      lr                    = lr,
      batch_size            = batch_size,
      n_aug_smp             = n_aug_smp,
      train_size            = train_size,
      seed                  = seed,
      use_gpu               = train_gpu > 0,
      docker                = docker,
      disk_size             = train_disk_size,
      mem_size              = train_mem_size,
      cpu                   = train_cpu,
      gpu                   = train_gpu
  }

  # ── 03a: Evaluate checkpoints, select model_order ─────────────────────────
  call Evaluate {
    input:
      anndata_h5ad           = preprocessed_h5ad,
      checkpoints_manifest   = TrainMixVAE.checkpoints_manifest,
      model_tar              = TrainMixVAE.model_tar,
      k_select_thr           = k_select_thr,
      batch_size             = batch_size,
      seed                   = seed,
      docker                 = docker,
      disk_size              = eval_disk_size,
      mem_size               = eval_mem_size,
      cpu                    = eval_cpu
  }

  output {
    # ── Human-review outputs ─────────────────────────────────────────────────
    # Download evaluation_results.json, inspect the consensus heatmap and
    # K-selection curve, then supply evaluation_results_json to MMIDAS_Analyze.
    File   evaluation_results_json = Evaluate.evaluation_results_json
    File   checkpoints_manifest    = TrainMixVAE.checkpoints_manifest
    File   model_tar               = TrainMixVAE.model_tar
    Array[File] evaluation_figures = Evaluate.evaluation_figures

    # ── Optional augmenter ───────────────────────────────────────────────────
    File?  augmenter_checkpoint    = TrainAugmenter.augmenter_checkpoint

    String pipeline_version_out    = pipeline_version
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task: TrainAugmenter  (02a — optional)
#
# Trains the UDAGAN VAE-GAN augmenter and writes:
#   RNA_augmenter_<tag>_<ts>.pth  — model checkpoint
#   loss_curve.png                — training loss
#   augmenter_manifest.json
# ──────────────────────────────────────────────────────────────────────────────
task TrainAugmenter {
  input {
    File   anndata_h5ad
    String tag
    Int    z_dim
    Int    noise_dim
    Int    fc_dim
    Int    n_epoch
    Int    batch_size
    Boolean use_gpu
    String docker
    Int    disk_size
    Int    mem_size
    Int    cpu
    Int    gpu
  }

  parameter_meta {
    anndata_h5ad: "Preprocessed AnnData .h5ad file (output of MMIDAS_DataPrep)."
    tag:          "Short label embedded in the checkpoint filename."
    z_dim:        "Latent space dimension of the augmenter encoder."
    noise_dim:    "Additive noise dimension."
    fc_dim:       "Hidden layer width."
    n_epoch:      "Number of training epochs."
    batch_size:   "Mini-batch size."
    use_gpu:      "Whether to pass --cuda to the script."
    docker:       "MMIDAS Docker image URI."
    disk_size:    "Boot disk size in GB."
    mem_size:     "Memory in GB."
    cpu:          "Number of CPU cores."
    gpu:          "Number of GPUs (0 = CPU only)."
  }

  command <<<
    set -euo pipefail
    mkdir -p out

    python3 /usr/local/02a_train_augmenter.py \
      --anndata_h5ad  "~{anndata_h5ad}" \
      --output_dir    out \
      --tag           "~{tag}" \
      --z_dim         ~{z_dim} \
      --noise_dim     ~{noise_dim} \
      --fc_dim        ~{fc_dim} \
      --n_epoch       ~{n_epoch} \
      --batch_size    ~{batch_size} \
      ~{if use_gpu then "--cuda" else ""}
  >>>

  output {
    File augmenter_checkpoint = glob("out/RNA_augmenter_*.pth")[0]
    File augmenter_manifest   = "out/augmenter_manifest.json"
    File loss_curve           = "out/loss_curve.png"
  }

  runtime {
    docker:       docker
    disks:        "local-disk ~{disk_size} HDD"
    memory:       "~{mem_size} GiB"
    cpu:          cpu
    gpuCount:     gpu
    gpuType:      "nvidia-tesla-t4"
    preemptible:  1
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task: TrainMixVAE  (02b)
#
# Trains the cpl-mixVAE with iterative pruning and writes:
#   model/cpl_mixVAE_model_before_pruning_<ts>.pth
#   model/cpl_mixVAE_model_after_pruning_<N>_<ts>.pth  (one per pruning round)
#   checkpoints_manifest.json
#
# The model/ directory is tarred so all checkpoints travel as a single File
# to the Evaluate task, which needs them all to build the summary dict.
# ──────────────────────────────────────────────────────────────────────────────
task TrainMixVAE {
  input {
    File    anndata_h5ad
    File?   augmenter_checkpoint
    Int     n_categories
    Int     state_dim
    Int     n_arm
    Int     latent_dim
    Int     fc_dim
    Float   x_drop
    Float   s_drop
    Float   temp
    Float   tau
    Float   beta
    Float   lam
    Float   lam_pc
    String  training_mode
    Int     n_epoch
    Int     n_epoch_p
    Float   min_con
    Int     max_prun_it
    Float   lr
    Int     batch_size
    Int     n_aug_smp
    Float   train_size
    Int     seed
    Boolean use_gpu
    String  docker
    Int     disk_size
    Int     mem_size
    Int     cpu
    Int     gpu
  }

  parameter_meta {
    anndata_h5ad:         "Preprocessed AnnData .h5ad file."
    augmenter_checkpoint: "Optional UDAGAN augmenter .pth (output of TrainAugmenter)."
    n_categories:         "Upper-bound number of cell-type categories."
    state_dim:            "Continuous state latent variable dimension."
    n_arm:                "Number of coupled encoder arms."
    latent_dim:           "Low-dimensional embedding dimension."
    fc_dim:               "Hidden layer width."
    x_drop:               "Input-layer dropout probability."
    s_drop:               "State-variable dropout probability."
    temp:                 "Gumbel-softmax temperature."
    tau:                  "Softmax temperature for categorical variable."
    beta:                 "KL divergence regularization coefficient."
    lam:                  "Coupling factor between arms."
    lam_pc:               "Coupling factor for the reference arm."
    training_mode:        "Reconstruction loss: MSE or ZINB."
    n_epoch:              "Training epochs before pruning."
    n_epoch_p:            "Training epochs per pruning iteration."
    min_con:              "Minimum inter-arm consensus to retain a category."
    max_prun_it:          "Maximum number of pruning iterations."
    lr:                   "Adam optimizer learning rate."
    batch_size:           "Mini-batch size."
    n_aug_smp:            "Number of augmented samples per real sample (0 = disabled)."
    train_size:           "Fraction of cells used for training (rest = validation)."
    seed:                 "Random seed."
    use_gpu:              "Whether to pass --cuda."
    docker:               "MMIDAS Docker image URI."
    disk_size:            "Boot disk size in GB."
    mem_size:             "Memory in GB."
    cpu:                  "Number of CPU cores."
    gpu:                  "Number of GPUs (0 = CPU only)."
  }

  command <<<
    set -euo pipefail
    mkdir -p out

    python3 /usr/local/02b_train_mixvae.py \
      --anndata_h5ad    "~{anndata_h5ad}" \
      --output_dir      out \
      ~{if defined(augmenter_checkpoint) then "--augmenter_checkpoint \"" + select_first([augmenter_checkpoint]) + "\"" else ""} \
      --n_categories    ~{n_categories} \
      --state_dim       ~{state_dim} \
      --n_arm           ~{n_arm} \
      --latent_dim      ~{latent_dim} \
      --fc_dim          ~{fc_dim} \
      --x_drop          ~{x_drop} \
      --s_drop          ~{s_drop} \
      --temp            ~{temp} \
      --tau             ~{tau} \
      --beta            ~{beta} \
      --lam             ~{lam} \
      --lam_pc          ~{lam_pc} \
      --training_mode   ~{training_mode} \
      --n_epoch         ~{n_epoch} \
      --n_epoch_p       ~{n_epoch_p} \
      --min_con         ~{min_con} \
      --max_prun_it     ~{max_prun_it} \
      --lr              ~{lr} \
      --batch_size      ~{batch_size} \
      --n_aug_smp       ~{n_aug_smp} \
      --train_size      ~{train_size} \
      --seed            ~{seed} \
      ~{if use_gpu then "--cuda" else ""}

    # Tar the model directory so all checkpoints travel as one File
    tar -czf model.tar.gz -C out model/
  >>>

  output {
    File checkpoints_manifest = "out/checkpoints_manifest.json"
    File model_tar            = "model.tar.gz"
  }

  runtime {
    docker:       docker
    disks:        "local-disk ~{disk_size} HDD"
    memory:       "~{mem_size} GiB"
    cpu:          cpu
    gpuCount:     gpu
    gpuType:      "nvidia-tesla-t4"
    preemptible:  0
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Task: Evaluate  (03a)
#
# Evaluates all checkpoints, runs K_selection to choose model_order, and writes:
#   evaluation_results.json  ← KEY HAND-OFF FILE for human review
#   consensus_T1_vs_T2_K_<K>.png
#   norm_consensus_T1_vs_T2_K_<K>.png
#   state_mu_K_<K>_arm_<arm>.png
#   summary_performance_K_<K>_narm_<N>.p
#
# Human review step: inspect evaluation_results.json and figures to confirm
# model_order is biologically sensible before submitting MMIDAS_Analyze.
# ──────────────────────────────────────────────────────────────────────────────
task Evaluate {
  input {
    File   anndata_h5ad
    File   checkpoints_manifest
    File   model_tar
    Float  k_select_thr
    Int    batch_size
    Int    seed
    String docker
    Int    disk_size
    Int    mem_size
    Int    cpu
  }

  parameter_meta {
    anndata_h5ad:         "Preprocessed AnnData .h5ad file."
    checkpoints_manifest: "checkpoints_manifest.json from TrainMixVAE."
    model_tar:            "Tarball of the model/ checkpoint directory from TrainMixVAE."
    k_select_thr:         "Minimum consensus threshold for K_selection(). Default: 0.95."
    batch_size:           "Batch size for inference."
    seed:                 "Random seed (must match 02b)."
    docker:               "MMIDAS Docker image URI."
    disk_size:            "Boot disk size in GB."
    mem_size:             "Memory in GB."
    cpu:                  "Number of CPU cores."
  }

  command <<<
    set -euo pipefail
    mkdir -p out

    # Restore the model checkpoints into the same relative path the manifest expects
    tar -xzf "~{model_tar}" -C out/

    # The manifest contains absolute paths from the TrainMixVAE task's working
    # directory. Rewrite them to point to our local out/ directory.
    MANIFEST_LOCAL="out/checkpoints_manifest.json"
    cp "~{checkpoints_manifest}" "$MANIFEST_LOCAL"
    TRAIN_OUTPUT_DIR=$(python3 -c "import json; print(json.load(open('$MANIFEST_LOCAL'))['output_dir'])")
    sed -i "s|${TRAIN_OUTPUT_DIR}|$(pwd)/out|g" "$MANIFEST_LOCAL"

    python3 /usr/local/03a_evaluate.py \
      --anndata_h5ad          "~{anndata_h5ad}" \
      --checkpoints_manifest  "$MANIFEST_LOCAL" \
      --output_dir            out \
      --k_select_thr          ~{k_select_thr} \
      --batch_size            ~{batch_size} \
      --seed                  ~{seed}
  >>>

  output {
    File        evaluation_results_json = "out/evaluation_results.json"
    Array[File] evaluation_figures      = glob("out/*.png")
  }

  runtime {
    docker:      docker
    disks:       "local-disk ~{disk_size} HDD"
    memory:      "~{mem_size} GiB"
    cpu:         cpu
    preemptible: 3
  }
}
