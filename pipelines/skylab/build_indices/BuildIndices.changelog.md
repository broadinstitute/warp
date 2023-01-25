# 2.0.1
2023-01-19 (Date of Last Commit)

* Added "Disk" to task runtime sections to support running on Azure

# 2.0.0

2022-12-20 (Date of Last Commit)

* Removed all tasks unrelated to Optimus
* Updated modify_gtf.py script to support references from Refseq and Gencode
* Updated the input files to accept gtf and fasta files directly
* 
# 1.0.1

2022-09-21 (Date of Last Commit)

* Docker image follows our guidelines
* Changed the type of biotypes from String to File so it localizes properly
* Changed the genome_fa to use the reference’s value instead of a modified_genome_fa that didn’t exist (which STAR was looking for and was then failing)

# 1.0.0

2022-02-01 (Date of Last Commit)

* Added modify_gtf.py and new docker for single nucleus smart-seq pipeline
* Update STAR to version 2.7.10a 
* Added biotypes as an input 

# 0.1.1

2021-11-15 (Date of Last Commit)

* Updated add-introns-to-gtf.py to use python3 instead of python2.

# 0.1.0

2021-05-03 (Date of Last Commit)

* Added a task to modify gtfs and fasta files and build indices for Single Nucleus Smart-seq pipeline
* Added a docker for the BuildStarSingleNucleus task.

# 0.0.1

2021-04-07 (Date of Last Commit)

* Added a star dockerfile to STAR version 2.7.8a, previously it was 2.5.3a


