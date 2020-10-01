# 2.1.1
2020-10-01

* Removed extra trailing slash in ouput directory from cloud to cloud copy job

# 2.1.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 2.0.2
2020-07-31

* Update various tasks to stop using phusion/baseimage:latest docker image (it has been removed).  Start using a Google-hosted base image in it's stead.

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
* Updated CopyFilesFromCloudToCloud to include the base_file_name in the contamination file.
* Fixed bug in CramToUnmappedBams.RevertSam where it was not reverting the OA tag
* Updated CramToUnmappedBams to properly use the output_map file to support testing.
* Renamed GermlineSingleSampleReferences to DNASeqSingleSampleReferences
* Updated shared tasks to support the new TargetedSomaticSingleSample pipeline


# 1.4
2020-03-05

* Added 'additional_disk' parameter to accommodate larger samples that have steps that run out of disk.

# 1.3
2019-12-03

* Modified the underlying ExomeGermlineSingleSample pipeline to use an up-to-date set of contamination resource files for VerifyBamId.  Further, the contamination resource files are subsetted by the target interval list.

# 1.2
Adjusted memory parameters to avoid Google's new e2 instances because there are not enough machines to satisfy our production use case.

# 1.1
This update is a result to an update to the ExomeGermlineSingleSample pipeline. 
We are jumping forward several versions of Picard, from version [2.18.27](https://github.com/broadinstitute/picard/releases/tag/2.18.27) to [2.20.4](https://github.com/broadinstitute/picard/releases/tag/2.20.4)

# 1.0
Initial release of the ExternalExomeReprocessing pipeline
