# 1.3.0
2022-05-24 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.4.1 to address log4j vulnerabilities

# 1.2.2
2022-03-11 (Date of Last Commit)

* Updated to Picard version 2.26.10 in CramToUnmappedBams subworkflow

# 1.2.1
2021-11-10

* Added Xmx flag (maximum heap size) to all tasks with java commands except mark duplicates

# 1.2.0
2021-07-19

* Changed GoTC image to Samtools specific image in CramToUnmappedBams and Utilities
* Changed GoTC image to GATK specific image in GermlineVariantDiscovery

# 1.1.0
2021-03-23

* Added a compression level argument to mark duplicates task and changed default disk and memory values

# 1.0.1
2021-03-02

* GDCWholeGenomeSomaticSingleSample pipeline has been moved out of the beta-pipelines folder and into the pipelines folder

# 1.0.0
2020-12-16

* First release of the GDCWholeGenomeSomaticSingleSample pipeline, pipeline produces identical results to GDC Alignment pipeline when run with example inputs