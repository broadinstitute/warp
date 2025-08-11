# 1.2.3
2025-08-11 (Date of Last Commit)

* Reorganized pipelines into the wdl directory

# 1.2.2
2024-11-04 (Date of Last Commit)

* Updated to GATK version 4.6.1.0

# 1.2.1
2024-09-10 (Date of Last Commit)

* Update to BGE filtering options in JointCalling; this has no effect on this pipeline

# 1.2.0
2024-09-06 (Date of Last Commit)

* Updated to GATK version 4.6.0.0
* Updated how no-calls are represented in the output VCFs (0/0 -> ./.) which also changes some annotations in the VCF

# 1.1.7
2023-12-18 (Date of Last Commit)

* Updated to GATK version 4.5.0.0.

# 1.1.6
2023-02-06 (Date of Last Commit)

* Updated VETS filtering pipeline to GATK version 4.5.0.0. Does not affect outputs.

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