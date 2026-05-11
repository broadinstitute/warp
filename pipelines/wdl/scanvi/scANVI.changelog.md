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

