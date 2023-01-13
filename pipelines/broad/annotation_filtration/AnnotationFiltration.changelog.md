# 1.2.4
2022-11-09 (Date of Last Commit)

* Updated to GATK version 4.3.0.0

# 1.2.3
2022-04-14 (Date of Last Commit)

* Updated to Picard version 2.26.10 and GATK version 4.2.6.1 to address log4j vulnerabilities

# 1.2.2
2022-01-12

* Fixed error in disk syntax in GatherFiltrationReport runtime block

# 1.2.1
2021-11-10

* Added Xmx flag (maximum heap size) to all tasks with java commands

# 1.2.0
2020-09-25

* Updated GATK docker image for all tasks to [GATK 4.1.8.0](https://github.com/broadinstitute/gatk/releases/tag/4.1.8.0).
    * Numerous bug fixes and improvements

# 1.1.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 1.0
Initial release of the AnnotationFiltration pipeline
