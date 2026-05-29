# aou_9.1.0
2026-05-20 (Date of Last Commit)

* Removed non-sharded option for step 3 (vcf merge) and unnecessary merge checks 

* Removed 'skip_summary' option from 'annotate_coverage' task

* Updated docker image references from `aou-mitochondrial-combine-vcfs-covdb` to `aou-mito-hail-processing:1.0.0` as part of mito docker consolidation

* Moved dockers and scripts from Warp to Warp-tools 

* Renamed inputs from "step3" to "vcf_merge" 

* Renamed workflow from "mt_coverage_merge" to "MitochondriaMerge"

* Updated per-position median calculation in `annotate_coverage` to match `hl.median` behavior: for even N, the median is now the floor of the average of the two middle values, rather than the lower of the two middle values

* Updated coverage HDF5 dtype from inferred uint16 to a fixed uint32 default to prevent silent truncation of high-coverage positions; added an explicit overflow guard that raises an error if any coverage value exceeds the dtype maximum

# aou_9.0.0
2026-03-23 (Date of Last Commit)

* Version of the pipeline use to process AoU v9 data


# aou_9_beta
2025-10-31 (Date of Last Commit)

* Added support for optional subsetting of inputs using a Terra data table TSV
* Adjusted resources for full size dataset


# beta release
2025-08-07 (Date of Last Commit)

* This is a working version of the mitochondria_single_sample and the starting point for the pipeline dev team