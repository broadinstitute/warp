# prototype
2026-09-28 (Date of Last Commit)

* Initial prototype of the phased-VCF-input counterpart to sc1_create_basis_mt_from_vds.py
* Imports phased per-chromosome VCF(s) directly via hl.import_vcf instead of reading a VDS -- no LGT/LA densification needed since VCF genotypes are already dense
* Annotates and filters against all 3 UCSC BED files (exome, clinvar, acaf_threshold), recomputes variant_qc/info (AC, AF, AN, homozygote_count)
* Test-only helpers not present in the VDS version: hl.filter_intervals-based --test_intervals, --sample_ids subsetting, enforced gs:// bucket paths, and no-clobber output checks
* Validated interactively against real phased chr21 data -- see sc1_create_basis_mt_from_vcf_validation_notes.md for the test region, methodology, and results
* Not yet wired into a WDL
