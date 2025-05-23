# aou-9.0.0
2025-05-20 (Date of Last Commit)
* First version of the RNAseq AoU pipeline
* Added pipeline_version to inputs
* Updated the STAR task docker from gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10 to gcr.io/broad-cga-francois-gtex/gtex_rnaseq@sha256:80c2db3cec3c08630237665e2d2f044f065022e0bbf7a62d0765f51f811818e2 to fix a STAR error for an unrecognized parameter
* Added the remove_ids task from the GTEx GitHub repository
