---
sidebar_position: 4
slug: /All_of_Us/Ancestry_Analysis/run_ancestry
title: Run Ancestry Inference
className: aou-doc-page
---

<div className="aou-folder-text">

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [aou_9.0.0](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_ancestry.changelog.md) | May, 2025 | WARP Pipelines | [File an issue](https://github.com/broadinstitute/warp/issues) |

## Introduction to the Run Ancestry Inference workflow

[`run_ancestry`](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_ancestry.wdl) is a WDL workflow that trains and applies a PCA-based ancestry classifier. It uses a filtered training callset (e.g., HGDP) to build PCA loadings and a Random Forest model, then projects target samples and predicts ancestry labels.

The workflow runs three major stages: training PCA model generation (`create_hw_pca_training`), ancestry prediction (`call_ancestry`), and visualization (`plot_ancestry`). Outputs include prediction tables, Hail Table tarballs, model artifacts, and interactive PCA plots.

## Quickstart table

| Pipeline Feature | Description | Source |
| :--: | :-- | :--: |
| Analysis type | PCA projection + Random Forest ancestry inference | |
| Workflow language | WDL 1.0 | [openWDL](https://github.com/openwdl/wdl) |
| Data input file format | VCF BGZF + Tabix index + optional metadata TSV | |
| Data output file format | TSV, `.tar.gz` Hail tables, `.pkl`, `.html` | |
| Primary software | Hail, scikit-learn, bokeh | [Hail](https://hail.is/), [scikit-learn](https://scikit-learn.org/), [Bokeh](https://bokeh.org/) |

## Set-up

### Run Ancestry Inference installation and requirements

The workflow code can be downloaded by cloning the [WARP GitHub repository](https://github.com/broadinstitute/warp). For the latest release, please see the [run_ancestry changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_ancestry.changelog.md).

The pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant workflow management system.

## Inputs

### Input descriptions

| Input variable name | Description | Type |
| --- | --- | --- |
| `hq_variants_intersection` | Sites-only VCF of variants shared between training and input datasets. | File |
| `hq_variants_intersection_idx` | Index for `hq_variants_intersection`. | File |
| `merged_vcf_shards` | Merged input-data VCF filtered to shared HQ sites. | File |
| `merged_vcf_shards_idx` | Index for `merged_vcf_shards`. | File |
| `filtered_training_set` | Training full VCF filtered to shared HQ sites. | File |
| `filtered_training_set_idx` | Index for `filtered_training_set`. | File |
| `hgdp_metadata_file_in` | *(Optional)* Training sample metadata TSV. Defaults to the public gnomAD HGDP+1KG metadata file. | File? |
| `final_output_prefix` | Prefix applied to all output artifacts. | String |
| `other_cutoff_in` | *(Optional)* Probability threshold for assigning ancestry label `oth` (other). Default: `0.75`. | Float? |
| `num_pcs` | Number of principal components used in training/projection. Default: `16`. | Int |

## Run Ancestry Inference tasks and tools

The workflow runs model training, prediction, and plotting tasks in sequence.

1. [Create training PCA artifacts](#1-create-training-pca-artifacts)
2. [Predict ancestry labels](#2-predict-ancestry-labels)
3. [Generate PCA plots](#3-generate-pca-plots)

To see specific tool parameters, select the task WDL link in the table; then view the `command {}` section of the task in the WDL script.

| Task name and WDL link | Tool | Software | Description |
| --- | --- | --- | --- |
| [create_hw_pca_training](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_ancestry.wdl) | Hail PCA | `hailgenetics/hail:0.2.67` | Runs HWE-normalized PCA on training data, writes loadings and labeled training PCA outputs. |
| [call_ancestry](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_ancestry.wdl) | Hail + RandomForestClassifier | `hailgenetics/hail:0.2.67` | Projects target samples, trains/applies RF classifier, and writes prediction + model outputs. |
| [plot_ancestry](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_ancestry.wdl) | Hail + bokeh | `hailgenetics/hail:0.2.67` | Produces interactive ancestry PCA HTML plots (raw and `oth`-adjusted labels). |

### 1. Create training PCA artifacts

Builds PCA features and loadings from the training set, joins population labels, and writes tarred Hail table artifacts plus eigenvalues.

### 2. Predict ancestry labels

Projects input samples into training PCA space, applies Random Forest classification, writes TSV predictions, Hail tables, and classifier pickle output.

### 3. Generate PCA plots

Uses prediction outputs to render interactive bokeh plots showing ancestry assignments on principal components.

## Outputs

| Output variable name | Filename, if applicable | Output format and description |
| --- | --- | --- |
| `results_tsv` | `<final_output_prefix>.ancestry_preds.tsv` | Tab-delimited prediction table containing sample IDs, predicted labels, probabilities, and PCA features. |
| `results_ht` | `<final_output_prefix>.ancestry_preds.ht.tar.gz` | Tarred Hail Table with projected PCA scores for input samples. |
| `results_loadings_ht` | `<final_output_prefix>_loadings.ht.tar.gz` | Tarred Hail Table containing PCA loadings used for projection. |
| `pred_plot` | `<final_output_prefix>.preds.html` | Interactive ancestry PCA plot using direct model predictions. |
| `pred_oth_plot` | `<final_output_prefix>.preds_oth.html` | Interactive PCA plot after applying `oth` thresholding. |
| `training_pca_labels_ht_tsv` | `<final_output_prefix>_training_pca.tsv` | Training PCA table with sample labels. |
| `eigenvalues_txt` | `<final_output_prefix>_eigenvalues.txt` | PCA eigenvalues text file used in model training. |
| `classifier_pkl` | `<final_output_prefix>_rf_classifier.pkl` | Pickled Random Forest classifier artifact. |

## Versioning

All `run_ancestry` releases are documented in the [changelog](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/run_ancestry.changelog.md).

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.

</div>
