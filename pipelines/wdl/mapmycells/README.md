# MapMyCells Pipeline

## Summary

MapMyCells is a cloud-optimized WDL workflow that performs highly accurate **cell type label transfer** for single-cell or single-nucleus RNA-seq (scRNA-seq/snRNA-seq) data. The pipeline is powered by the Allen Institute's `cell_type_mapper` library, which relies on an extremely fast hierarchical nearest-centroid algorithm to map query data against massive precomputed reference atlases (e.g., the 4M-cell Allen Brain Cell Atlas).

Unlike iterative machine learning methods (such as scANVI or Seurat), MapMyCells maps queries to a statistically precomputed representation of the reference taxonomy. This gives it two distinct advantages:
1. It scales instantly and maps hundreds of thousands of cells in minutes on a standard CPU.
2. It flawlessly resolves exceptionally rare subpopulations and extremely granular hierarchical taxonomies.

## Reference Atlases

The pipeline supports three modes of mapping controlled by the `reference_atlas` input:
* `"Human_MTG"` (Default) — Maps the query against the Human Middle Temporal Gyrus (MTG) taxonomy.
* `"Mouse_WMB"` — Maps the query against the massive Whole Mouse Brain (WMB) / Allen Brain Cell (ABC) taxonomy.
* `"Custom"` — Allows the user to provide their own `.h5` reference statistics and `.json` markers.

## Inputs

The MapMyCells pipeline takes an unannotated `query_h5ad` file and outputs `.csv` and `.json` predictions.

| Parameter | Type | Description |
|-----------|------|-------------|
| `MapMyCells.query_h5ad` | File | **Required.** The input AnnData file containing raw scRNA-seq counts. |
| `MapMyCells.input_id` | String | **Required.** An identifier prepended to all output filenames. |
| `MapMyCells.reference_atlas` | String | Defines which precomputed atlas is used. Options: `"Human_MTG"`, `"Mouse_WMB"`, `"Custom"`. Default: `"Human_MTG"`. |
| `MapMyCells.custom_precomputed_stats` | File? | **Required if Custom.** An HDF5 file containing the precomputed hierarchical taxonomy stats. |
| `MapMyCells.custom_query_markers` | File? | *Optional.* JSON file containing predefined marker genes. Only used when `reference_atlas` is `"Custom"` — Human_MTG and Mouse_WMB use their own baked-in markers (or none, for Human_MTG). |
| `MapMyCells.custom_gene_mapping_db` | File? | *Optional.* SQLite DB for translating gene symbols to Ensembl IDs. Defaults to Allen Institute's pre-built db (`gs://broad-gotc-test-storage/mapmycells/mmc_gene_mapper.2025-08-04.db`, ~15GB) for every atlas; override only if you need a different mapping. |
| `MapMyCells.algorithm` | String | Type assignment algorithm. Options: `"hierarchical"`, `"hann"`. Default: `"hierarchical"`. |

## Outputs

The task generates two output files:
* `~{input_id}_mapmycells_results.csv`: Contains the hierarchical mapping assignment for each cell (class, subclass, supertype, cluster) along with bootstrapping probabilities.
* `~{input_id}_mapmycells_extended.json`: A verbose JSON output detailing the tree traversal, probabilities, and mapped pathways for each individual cell.

## Usage

You can find example inputs in the `test_inputs/` folder:
* `test_inputs/Scientific/test_human.json` — Maps using `"Human_MTG"`
* `test_inputs/Scientific/test_mouse.json` — Maps using `"Mouse_WMB"`

The Human_MTG and Mouse_WMB reference-atlas assets are baked into the pipeline's docker image (see [warp-tools/3rd-party-tools/mapmycells](https://github.com/broadinstitute/warp-tools/tree/develop/3rd-party-tools/mapmycells)), so you do not need to download the reference atlases or pass them into the workflow yourself. Simply specify the `reference_atlas` and provide your `query_h5ad`!

## Citing

This pipeline wraps the Allen Institute's `cell_type_mapper` library. If you use it in your research, please cite:

* Daniel SF, Lee C, Mollenkopf T, et al. *High-performance mapping of unlabeled cell-by-gene data to reference brain taxonomies.* bioRxiv (2026). doi: [10.64898/2026.03.06.710160](https://doi.org/10.64898/2026.03.06.710160) ([PMID: 41958981](https://pubmed.ncbi.nlm.nih.gov/41958981/))
* Source code: [github.com/AllenInstitute/cell_type_mapper](https://github.com/AllenInstitute/cell_type_mapper)
