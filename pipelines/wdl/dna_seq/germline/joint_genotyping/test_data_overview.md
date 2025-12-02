# Joint Genotyping Test Data
The following lists the data sets used for workflow validation, including workflow engineering (plumbing) tests and scientific validation tests.

* WARP tests (scientific and plumbing) run 5 configurations: Gnarly and GenotypeGVcFs with combinations of {gather, don't gather VCFs} and {serial, parallel VQSR}.

* All input GVCFs should be \"reblocked\" with latest version of ReblockGVCF, as in the single-sample pipeline.

### Exome Joint Calling
* Input GVCFs include 60 1000G samples plus SynDip, CEU trio, and NA19238 Yoruba trio daughter

### Whole Genome Germline (WGS) Joint Calling
* SynDip (WGS2) NA12878, husband (NA12877), parents YRI trio 60 more 1000G samples