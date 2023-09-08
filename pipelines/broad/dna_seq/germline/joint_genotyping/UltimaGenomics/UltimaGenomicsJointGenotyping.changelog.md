# 1.1.5
2023-09-08 (Date of Last Commit)

* Added option to hard filter sites outside of provided interval list to HardFilterAndMakeSitesOnlyVcf task

# 1.1.4
2023-06-29 (Date of Last Commit)

* Added extra_args input to the SplitIntervalList task to support the JointGenotypingTasks.wdl

# 1.1.3
2023-05-23 (Date of Last Commit)

* Made disk and memory available as inputs to the JointGenotypingTasks.wdl.

# 1.1.2
2023-03-01 (Date of Last Commit)

* Updated version of Ultima's variant calling package which is used to choose a threshold for filtering by mesauring F1 against a truth set
* Update includes a bug fix in concordance measurement that incorrectly discounted variants with wrong allele called

# 1.1.1
2022-11-29 (Date of Last Commit)

* Updated version of GATK to 4.3.0.0.
* Added additional boot disk to 4 tasks in UltimaGenomicsGermlineFilteringThreshold.wdl.

# 1.1.0
2022-08-24 (Date of Last Commit)

* Added filtering to UltimaGenomicsJointGenotyping pipeline. Output VCFs are now filtered.

# 1.0.0
2022-05-09 (Date of Last Commit)

* Initial release of the UltimaGenomicsJointGenotyping pipeline