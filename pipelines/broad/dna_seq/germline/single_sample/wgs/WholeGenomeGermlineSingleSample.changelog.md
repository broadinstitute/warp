# 2.3.3
2021-03-17

* Promoted VariantCalling to be a top-level workflow

# 2.3.2
2021-02-22

* Added SORTING_COLLECTION_SIZE_RATIO as an optional task input to MarkDuplicates

# 2.3.1
2020-12-21

* Passed an input bam index to several subworkflows, so the pipeline passes on singularity for sharded BQSR

# 2.3.0
2020-12-16

* Fixed error in relative import statement in Alignment subworkflow.
* Fixed syntax bug in Alignment task SamToFastqAndBwaMemAndMba

# 2.2.0
2020-10-20

* Updated GATK docker image for all tasks to [GATK 4.1.8.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.8.0).
    * Numerous bug fixes and improvements
* Updated Picard docker image for all tasks to [2.23.8](https://github.com/broadinstitute/picard/releases/tag/2.23.8).
* Updated samtools to version [1.11](https://github.com/samtools/samtools/releases/tag/1.11).  Primarily for improved compression of cram files.

# 2.1.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 2.0.1
2020-07-15

* Remove GetBWAVersion as a task and moved it to SamToFastqAndBwaMemAndMba

# 2.0 
2020-05-13

### Breaking changes to the structure of pipeline inputs. 
* Changes to the inputs included with the dna seq single sample references struct:
    * Removed 'fingerprint_genotypes_file' and 'fingerprint_genotypes_index' from bundle and made these optional pipeline inputs
    * Removed 'haplotype_scatter_count' and 'break_bands_at_multiples_of' from bundle and added these to a separate 'VariantCallingScatterSettings' struct
    * Added 'haplotype_database_file' to the references bundle as a non-optional file
#### Additional changes
* Renamed GermlineSingleSampleReferences to DNASeqSingleSampleReferences
* Updated shared tasks to support the new TargetedSomaticSingleSample pipeline

# 1.4
2020-03-05

* Added 'additional_disk' parameter to accommodate larger samples that have steps that run out of disk.

# 1.3
2019-12-03

* Modified the WholeGenomeGermlineSingleSample pipeline to use an up-to-date set of contamination resource files for VerifyBamId.
* Removed unused import statements from WholeGenomeGermlineSingleSample.wdl

# 1.2
Adjusted memory parameters to avoid Google's new e2 instances because there are not enough machines to satisfy our production use case.

# 1.1
This is a major update to the WholeGenomeGermlineSingleSample pipeline. We are jumping forward several versions of Picard, from version [2.16.0](https://github.com/broadinstitute/picard/releases/tag/2.16.0) to [2.20.4](https://github.com/broadinstitute/picard/releases/tag/2.20.4)
## Changes to Expect
For WGS we have observed the following changes in out test data:

### Crams
[OA](https://github.com/broadinstitute/picard/commit/fbb06096) tags have been added to Crams, PA tags are changed to OA.
Slight changes in the selection of a non-duplicate candidate from a duplicate set of reads.

### Metrics
- Alignment Summary Metrics and Readgroup Alignment Summary Metrics
  - Minor changes in the values of the `PCT_ADAPTER` metric.
- Bait Bias Detail Metrics
  - Minor changes in the values of the `FWD_CXT_REF_BASES`, `FWD_CXT_ALT_BASES`, `REV_CXT_REF_BASES`, and `REV_CXT_ALT_BASES` metrics.
- Bait Bias Summary Metrics
  - Minor changes in the values of the `WORST_POST_CXT_QSCORE` metric.
- Duplicate Metrics
  - Minor changes in the values of `READ_PAIR_OPTICAL_DUPLICATES` and `ESTIMATED_LIBRARY_SIZE`.
  - New metrics and additional bins included in the [histogram](https://github.com/broadinstitute/picard/commit/243a2ca4).
- Insert Size Metrics and Readgroup Insert Size Metrics
  - New metrics, `MODE_INSERT_SIZE` and `WIDTH_OF_95_PERCENT`, are now delivered. 
  - Small changes in values of the `MEAN_INSERT_SIZE`, `STANDARD_DEVIATION`, and `READ_PAIRS` metrics.
  - Differences in values in the histograms.
- Pre Adapter Detail Metrics
  - Small changes in values of the `PRO_REF_BASES`, `PRO_ALT_BASES`, `CON_REF_BASES`, `CON_ALT_BASES`, and `ERROR_RATE` metrics.
- Raw WGS metrics
  - New metric, `PCT_EXC_ADAPTER`, is now delivered. 
  - Small changes in values of the `MEAN_COVERAGE`, `PCT_EXC_BASEQ`, `PCT_EXC_DUPE`, `PCT_EXC_TOTAL`, `PCT_EXC_UNPAIRED`, and `SD_COVERAGE` metrics.
  - Differences in values in the histograms.
- WGS metrics
  - New metric, `PCT_EXC_ADAPTER`, is now delivered. 
  - Small changes in values of the `PCT_EXC_MAPQ` metric.
- Variant Calling Detail Metrics
  - Minor changes in the values of the `HET_HOMVAR_RATIO`, `TOTAL_HET_DEPTH`, `TOTAL_SNPS`, `NUM_IN_DB_SNP`, `NOVEL_SNPS`, `PCT_DBSNP`, `NOVEL_TITV`, `TOTAL_INDELS`, `NOVEL_INDELS`, `PCT_DBSNP_INDELS`, `NUM_IN_DB_SNP_INDELS`, `DBSNP_INS_DEL_RATIO`, `NOVEL_INS_DEL_RATIO`, `TOTAL_MULTIALLELIC_SNPS`, `TOTAL_COMPLEX_INDELS`, `NUM_IN_DB_SNP_COMPLEX_INDELS`, `SNP_REFERENCE_BIAS`, and `NUM_SINGLETONS` metrics.
- Variant Calling Summary Metrics
  - Minor changes in the values of the `TOTAL_SNPS`, `NUM_IN_DB_SNP`, `NOVEL_SNPS`, `PCT_DBSNP`, `NOVEL_TITV`, `TOTAL_INDELS`, `NOVEL_INDELS`, `PCT_DBSNP_INDELS`, `NUM_IN_DB_SNP_INDELS`, `DBSNP_INS_DEL_RATIO`, `NOVEL_INS_DEL_RATIO`, `TOTAL_MULTIALLELIC_SNPS`, `TOTAL_COMPLEX_INDELS`, `NUM_IN_DB_SNP_COMPLEX_INDELS`, `SNP_REFERENCE_BIAS`, and `NUM_SINGLETONS` metrics.

# 1.0
Initial release of the WholeGenomeGermlineSingleSample pipeline
