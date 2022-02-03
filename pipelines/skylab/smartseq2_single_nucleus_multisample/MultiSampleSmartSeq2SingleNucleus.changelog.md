# 1.2.0
2022-01-07 (Date of Last Commit)

* Fixed missing metadata issue in the loom file
# 1.1.2
2021-11-19 (Date of Last Commit)

* Updated STARsoloFastq to use 'output_bam_basename' to name the aligned bam. This is consistent with versions 4.2.7 and older. This change has no impact MultiSampleSmartSeq2SingleNucleus
# 1.1.1
2021-11-15 (Date of Last Commit)
* Updated remove-reads-on-junctions.py in the FeatureCounts.wdl to use python3 instead of python2.

# 1.1.0
2021-09-16 (Date of Last Commit)

* Removed the Smart-seq2 Single Nucleus workflow (SmartSeq2SingleNucleus.wdl) and changed the workflow tasks to run multiple samples in the same VM. This change is expected to make the pipeline faster and cheaper.
* Renamed the StarAlignFastq.StarAlignFastqPairedEnd task to StarAlign.StarAlignFastqMultisample

# 1.0.4
2021-09-10 (Date of Last Commit)

Added the option "--soloBarcodeReadLength 0" STARsoloFastq task to support alignment in Optimus. This change has no impact on MultiSampleSmartSeq2SingleNucleus

# 1.0.3
2021-09-02 (Date of Last Commit)

* Added a new StarSolo task for Optimus in the StarAlign.wdl. However, the same wdl
  contains other Star tasks that are used in the smartseq2 single nuclei for paired and 
  single stranged fastq files. As a result, the smartseq2 processing is not expected to 
  change. 

# 1.0.2
2021-08-02 (Date of Last Commit)

* Increased the version number to make new release tag for Dockstore 

# 1.0.1
2021-07-19 (Date of Last Commit)

* Updated SmartSeq2 to accommodate spaces in input_name

# 1.0.0

2021-05-17 (Date of First Commit)

* This is the first release of the Smart-seq2 Multi-Sample Single Nuclei workflow
