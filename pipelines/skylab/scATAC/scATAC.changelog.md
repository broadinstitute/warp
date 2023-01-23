# 1.3.1
2023-01-19 (Date of Last Commit)

* Added "Disk" to task runtime sections to support running on Azure

# 1.3.0
2022-09-23 (Date of Last Commit)

* Added disk, memory and cpu as task inputs. Added pipeline version as a string output.

# 1.2.4
2022-08-23 (Date of Last Commit)

* Remove an unused script in pytools docker image.

# 1.2.3
2022-08-18 (Date of Last Commit)

* Update AlignPairedEnd, SnapPre, SnapCellByBin tasks to use rebuilt snaptools docker image.

# 1.2.2
2022-08-16 (Date of Last Commit)

* Update MakeCompliantBAM and BreakoutSnap tasks to use a consolidated python utilities docker image.

# 1.2.1

2021-11-15 (Date of Last Commit)

* Updated breakoutSnap.py to use python3 instead of python2.

# 1.2.0

2020-12-23 (Date of Last Commit)

* Added input_id as an input to the scATAC workflow and used that to create unique output names

# 1.1.0

2020-08-24 (Date of Last Commit)

* Added bin_size_list as an optional input to the scATAC workflow with a default of 10,000 bp

# 1.0.1

2020-05-28 (Date of Last Commit)

* Changed the pipeline name from "snap-atac" to "scATAC"
* Updated the documentation

# 1.0

2019-11-08 (Date of Last Commit)

* Initial release of the snap-atac pipeline 
* The pipeline can process single-cell ATAC-seq fastq files that have been pre-tagged with cellular barcodes and produce an output snap file that includes peak summarization


