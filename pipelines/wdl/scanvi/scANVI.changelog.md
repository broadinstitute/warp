# 3.0.0
2026-10-29 (Date of Last Commit)

### Breaking changes

* Removed the placeholder sample-metadata columns that the scvi-scanvi container stamped onto every cell of the SCANVI_predictions output (biosample_id, donor_id, species, species__ontology_label, disease, disease__ontology_label, organ, organ__ontology_label, library_preparation_protocol, library_preparation_protocol__ontology_label, data_modality, sex). They were constants (sample_001, donor_001, Homo sapiens, brain, 10x multiome, ...) whatever the input, so they were wrong for most runs. The pipeline adds no sample metadata; sample metadata is joined downstream from the data archive's sample manifests
* Changed the sample inputs to arrays so one run labels one or more samples. input_ids and gex_h5ads, paired by position, replace input_id and gex_h5ad, and the optional atac_h5ads (one ATAC h5ad per GEX h5ad, or empty for GEX-only mode) replaces atac_h5ad. A single sample is a one-element array: an inputs JSON needs the brackets, while Terra wraps a single data-table value into an array
* Changed the outputs to arrays with one entry per sample, in input order. gex_annotated_h5ads, atac_annotated_h5ads (empty in GEX-only mode) and scanvi_models_out replace gex_annotated_h5ad, atac_annotated_h5ad and scanvi_model_out
* Changed the output file names to follow the input files: each input <name>.h5ad gives <name>_annotated.h5ad (replacing <input_id>_gex_annotated_matrix.h5ad and <input_id>_atac_annotated_matrix.h5ad), and a training run's models are saved as <gex name>_scanvi_model.tar.gz. The input files of a run must have distinct names, which the input check enforces; input_ids only records each sample's id in its outputs
* Removed the SCANVI predictions output (scanvi_predictions_h5ad); its labels and latent are in the annotated matrices. scanvi_models_out is now written only by a training run and is empty when a model is supplied
* Added label levels: a sample is labelled at one or more levels of the reference, such as class and subclass, each with its own SCANVI model trained from one shared SCVI model. ref_label_columns (a list, default the reference's subclass or final_annotation) replaces ref_label_column. Each level's labels go in obs['<level>'] (replacing obs['final_annotation']), its latent in obsm['X_scANVI_<level>'] and its UMAP in obsm['X_umap_<level>']; obsm['X_umap'] is a copy of the first level's UMAP, and uns['scanvi_levels'] lists the levels. With output_max_probability, each level's confidence goes in obs['<level>_max_probability'] (replacing obs['max_probability'])
* Changed the model archive (scanvi_model input and scanvi_models_out) to hold every model trained from one reference: reference.txt, levels.txt, the scVI model in scvi/ and one SCANVI model per level in scanvi/<level>/. An archive in the old single-model form is refused with a message naming the expected layout
* Removed bucket mode (input_bucket, gex_filename, atac_filename and ref_filename). Bucket mode could not run in the pinned image, which has no gsutil; list the files in gex_h5ads and atac_h5ads instead
* Changed the annotated GEX and ATAC matrices to the form Single Cell Portal reads: the counts move to .raw (float32), .X holds log2(CPM + 1) (float32; units are bits, the Allen Institute's convention, with the base recorded in uns), obsm gains each level's 30-dimensional scANVI latent and a 2-column UMAP fitted on the sample's own latent, and the duplicate counts layer and the placeholder final_annotation = Unknown column are dropped. Readers that took counts from .X or layers['counts'] must read .raw instead
* Removed the PreprocessFilter task and its nested inputs. Each sample's preprocessing now runs inside the label-transfer task, in the container's preprocess_query.py, with the same steps
* Changed the label-transfer task to run every sample of a run on one machine, one after another, instead of one machine per sample. Its default size drops from 32 CPUs and 120 GiB to 16 CPUs and 104 GiB, which Terra finds more readily; memory only has to fit the largest sample

### Additional changes

* Updated the pinned scvi-scanvi docker image to image 1.5.0, with scvi-tools 1.5.1, snapatac2 2.10.0, anndata 0.13.4, scipy 1.18.1 and torch 2.15.0 (previously scvi-tools 1.2.0, snapatac2 2.7.0, and scipy held below 1.14). SCVI/SCANVI training and prediction now run on the newer libraries, which may produce slightly different (not qualitatively different) label assignments and embeddings. Each output records the image's key library versions in uns['scanvi_versions']
* Added the optional seed input, which makes SCVI/SCANVI training reproducible; unset (the default) keeps training nondeterministic
* Updated model loading to extract a supplied scanvi_model archive with Python's safe data filter
* Added query mapping for models trained on the reference alone, such as those from warp-scanvi-training. When a supplied model does not know a sample's batch label, the sample is mapped onto it with scArches (load_query_data plus a short training on the sample alone), set by the new optional query_max_epochs input (default 100). A model this pipeline trained still loads as is
* Made ref_h5ad optional: it is needed only to train. A run with a supplied scanvi_model neither prepares nor loads the reference, so its cost no longer grows with the reference's size; the run fails at once when neither is given
* Added the PreprocessReference task, which prepares the reference once per training run and shares it with every sample, on a 16-CPU, 104 GiB machine with disk sized from the reference
* Added run provenance to every output h5ad: uns records scanvi_input_id (the sample's id), scanvi_input_file (the GEX file's name), scanvi_model (the model's name, which is the reference it was trained on: the reference's title, or its file name, plus the first 8 hex digits of the reference file's sha256) and scanvi_levels. A training run writes that name into each model archive as reference.txt, and a supplied model's reference.txt is carried forward. uns also records scanvi_versions, the key library versions of the image that wrote the file. The scanvi_ prefix keeps these apart from the input_id and pipeline_version an Optimus h5ad already carries
* Added a check that fails the run before any sample is processed when input_ids, gex_h5ads and atac_h5ads do not line up, and a check in the label-transfer task that every input file opens as an h5ad, the ids are unique and no two input files share a name, reporting every problem at once
* Changed disk_size to an optional input: when unset, the label-transfer disk is sized from the inputs (4 times the samples' input size, plus 3 times the reference size when training, plus the model and 50 GB)
* Added gzip compression to every output h5ad

# 2.1.0
2026-07-09 (Date of Last Commit)

* Added an optional batch_size input (default 128) that sets the SCVI/SCANVI minibatch size in both multiome and GEX-only modes. Lower it to fit a high-cardinality reference on a small GPU (SCANVI activation memory scales with batch_size x number of labels); raise it on a large-VRAM cloud GPU. Default 128 is scvi-tools' own default, so existing outputs are unchanged.
* Output metadata: added a `data_modality` obs field ("Simultaneous profiling of gene expression and open chromatin from the same cell.") and corrected the library-prep ontology label to "10x multiome" (`library_preparation_protocol` = "EFO_0030059", `library_preparation_protocol__ontology_label` = "10x multiome"), in the container's finalize_output.
* Reworked the label-transfer step: moved its logic into the container (new label_transfer_from_preprocessed.py entry point) so the WDL task is a one-line command, and split it into a GPU task (MultiomeLabelTransfer) and a CPU-only task (MultiomeLabelTransferCpu) selected by gpu_count. This makes gpu_count = 0 a real CPU-only run (e.g. supplied-model prediction) — Cromwell rejects a gpuCount of 0 and cannot conditionally omit the GPU runtime attributes from a single task, so the CPU path needs its own task. Both tasks run the identical command; scvi-tools auto-detects the accelerator. GPU runtime attributes also switched to portable camelCase (gpuType/gpuCount/nvidiaDriverVersion).
* Updated the pinned scvi-scanvi docker image to the sc-scanvi-gex-only build that includes the label_transfer_from_preprocessed.py entry point (and the batch_size / metadata changes above).

# 2.0.0
2026-06-30 (Date of Last Commit)

* MAJOR: the pipeline now emits a reusable model. Every run outputs its trained (or loaded) SCANVI model as scanvi_model_out — a saved model directory (.tar.gz, ~tens of MB, no bundled adata), in the same format the pipeline accepts as input. Reusing that model to run inference on similar data means reprocessing through the pipeline, so this is a qualitative change to the outputs (major per VersionAndReleasePipelines.md).
* Added an optional scanvi_model input: when provided, the label-transfer task loads it and predicts instead of training SCVI/SCANVI — valid when the model matches the incoming data/reference. Auto-detected (no flag).
* Exposed gpu_count (default 2) plus mem_size/nthreads/disk_size as inputs so a supplied-model prediction run can be right-sized to a small CPU box (gpu_count=0); default training/inference resources are unchanged.
* Added an optional output_max_probability input (default false). When true, every output h5ad gains a max_probability obs column holding the per-cell maximum SCANVI posterior probability (the confidence of the assigned label), computed from lvae.predict(soft=True). Default false leaves existing outputs unchanged.
* Added support for AIT (Allen Institute Taxonomy) schema reference atlases. AIT references are auto-detected (uns['schema_version'] + uns['hierarchy']) and adapted in PreprocessFilter: counts are materialized from .raw (AIT files have no .X) using the gene symbols from .var, the cell-type label is taken from a chosen taxonomy level, and the batch from a chosen obs column.
* Added optional ref_label_column and ref_batch_column inputs to select the reference label/batch columns. Defaults: AIT references use subclass/donor_id; PBMC-style references use final_annotation/batch (existing behavior unchanged).
* Added a genome input (hg38 default, mm10, mm39) for the ATAC cell-by-bin -> gene-activity conversion, enabling mouse multiome references; default hg38 keeps existing human runs unchanged.

# 1.1.0
2026-06-27 (Date of Last Commit)

* Made the ATAC h5ad input optional. When no ATAC h5ad is provided (direct-file mode) or no ATAC object is found in the input bucket (bucket mode), the pipeline auto-detects GEX-only mode and trains/annotates from the reference atlas using gene expression and reference data alone, without using ATAC.
* Added run_gex_only_model to the scvi-scanvi container to train SCVI/SCANVI on GEX + reference only
* PreprocessFilter skips ATAC loading, barcode reindexing, shared-barcode subsetting, and gene-activity conversion in GEX-only mode
* Made preprocessed_atac_activity_h5ad and atac_annotated_h5ad optional outputs (null in GEX-only mode)
* Added a GEX-only example inputs file (example_inputs/scANVI.gex_only.json)
* Added an optional max_epochs input that caps SCVI/SCANVI training epochs in both multiome and GEX-only modes (defaults to the container value of 500 when unset)

# 1.0.0
2026-04-17 (Date of Last Commit)

* Initial release of the ScviScanvi pipeline for cell type label transfer on Multiome data
* Integrated SCVI and SCANVI models using the scvi-scanvi docker image (1.0.0-1.2-1756234975)
* Test GPU support for accelerated model training
* Outputs SCANVI predictions, annotated ATAC, and annotated GEX h5ad files
* Added PreprocessFilter task to handle all h5ad preprocessing and filtering as a separate CPU-only step before GPU model training
* Moved column patching (star_IsCell, gex_barcodes), GEX cell/gene filtering, ATAC barcode reindexing, shared barcode subsetting, batch/modality tagging, and gene activity matrix conversion out of the monolithic MultiomeLabelTransfer task
* Simplified MultiomeLabelTransfer task to accept preprocessed h5ad files directly; removed inline Python patching and bucket resolution logic
* Hardcode GPU configuration in runtime block for MultiomeLabelTransfer; no GPU needed for PreprocessFilter
* added test_cuda.wdl to the tasks directory

