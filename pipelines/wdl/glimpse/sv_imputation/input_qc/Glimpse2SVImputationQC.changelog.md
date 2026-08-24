# 0.0.1
2026-08-24 (Date of Last Commit)

* Initial release of pipeline to perform QC checks on inputs to the SV Imputation pipeline using GLIMPSE2.
* Checks include:
  - Input manifest:
    - file accessibility/existence (does not allow RP bucket)
    - presence of both gvcf_path and gvcf_index_path columns/column headers
    - each row has matching basename between the gvcf and its index
    - no duplicate paths in different rows
    - gvcf files have '.vcf.gz' or '.gvcf.gz' extension and index files have 'vcf.gz.tbi' extension
  - Individual GVCF files:
    - no input file exceeds a defined cutoff size (default 10GB)
    - contigs must match ref dict, including length values
    - header must contain PL and GT FORMAT annotations
    - header must contain multiple GVCFBlock lines
    - must be single sample
  - Uniqueness of sample identifiers extracted from GVCFs
