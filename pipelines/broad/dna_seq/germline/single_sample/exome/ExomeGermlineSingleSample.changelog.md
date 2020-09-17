# 2.1.0
2020-08-18

* Added 'allowNestedInputs: true' metadata parameter to wdl to support Cromwell version 52

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

* Modified the ExomeGermlineSingleSample pipeline to use an up-to-date set of contamination resource files for VerifyBamId.  Further, these contamination resource files are subsetted by the target interval list.
* Removed unused import statements from ExomeGermlineSingleSample.wdl

# 1.2
Adjusted memory parameters to avoid Google's new e2 instances because there are not enough machines to satisfy our production use case.

# 1.1
This is an update to the ExomeGermlineSingleSample pipeline. We are jumping forward several versions of Picard, from version [2.18.27](https://github.com/broadinstitute/picard/releases/tag/2.18.27) to [2.20.4](https://github.com/broadinstitute/picard/releases/tag/2.20.4)
## Changes to Expect
For Exomes we have observed the following changes in out test data:

### Crams
[OA](https://github.com/broadinstitute/picard/commit/fbb06096) tags have been added to Crams, PA tags are changed to OA.

### Metrics
- Alignment Summary Metrics and Readgroup Alignment Summary Metrics
  - Minor changes in the values of the `PCT_ADAPTER` metric.
- Duplicate Metrics
  - Minor changes in the values of `READ_PAIR_OPTICAL_DUPLICATES `and `ESTIMATED_LIBRARY_SIZE`.
  - New values and additional bins included in the histogram.
- Hybrid Selection Metrics
  - New metrics, `MIN_TARGET_COVERAGE`, `PCT_EXC_ADAPTER`, and `PF_BASES`, are now delivered. 
  - Small changes in values of the `ZERO_CVG_TARGETS_PCT ` metric.

# 1.0
Initial release of the ExomeGermlineSingleSample pipeline
