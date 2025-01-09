# 1.1.3
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 1.1.2
2022-04-14 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.1.1
2021-11-10

* Added Xmx flag (maximum heap size) to all tasks with java commands

# 1.1.0
2021-07-19

* Changed GoTC image to Samtools specific image in CramToUnmappedBams and Utilities
* Changed GoTC image to GATK specific image in GermlineVariantDiscovery

# 1.0.1
2021-02-08

* Calculate java memory value from the optional memory input value for CramToUnmappedBams java tasks

# 1.0.0
2021-02-02

### Initial release of the CramToUnmappedBams pipeline
* Pipeline to create unmapped bams from an aligned cram or bam file
