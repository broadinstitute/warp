# 5.9.2
2024-11-15 (Date of Last Commit)

* Added bam validation in the StarSoloFastq task; this does not affect the outputs of the pipeline

# 5.9.1
2024-11-12 (Date of Last Commit)

* Renamed the ATAC workflow library metric percent_target to atac_percent_target for compatibility with downstream tools
* Added more disk and memory to the JoinBarcodes task

# 5.9.0
2024-10-21 (Date of Last Commit)

* Updated the tabix flag in JoinMultiomeBarcodes task in H5adUtils.wdl to use CSI instead of TBI indexing, which supports chromosomes larger than 512 Mbp; this task changes the format for the ATAC fragment file index 
* Renamed the fragment file index from atac_fragment_tsv_tbi to atac_fragment_tsv_index


# 5.8.0
2024-10-23 (Date of Last Commit)

* Updated the workflow to include a new expected_cells input parameter describing the number of cells used as input to the library preparation; this is passed to both the ATAC workflows and Optimus workflows and the default is set to 3000 cells
* Updated the ATAC library CSV and the Gene Expression library CSV to be consistent in file naming convention and to have similar case for metric names
* Added a new metric to the ATAC library CSV to calculate percent_target, which is the number of estimated cells by SnapATAC2 divided by expected_cells input
* Updated the ATAC workflow so that the output fragment file is bgzipped by default
* Updated memory settings for PairedTag; does not impact the Multiome workflow

# 5.7.1
2024-10-18 (Date of Last Commit)

* Removed the underscore of the NHashID in the ATAC library metrics CSV to match the gene expression library metrics

# 5.7.0
2024-09-24 (Date of Last Commit)
* Added a python implementation of DoubletFinder to calculate doublet scores in gene expression data; percent doublets are now available as a library-level metric and individual doublet scores for cell barcodes are in the h5ad
* Updated gene_names in the final h5ad to be unique

# 5.6.1
2024-09-11 (Date of Last Commit)
* Updated warp-tools docker which added create_h5ad_snss2.py to the docker image. This change does not affect the Multiome pipeline

# 5.6.0
2024-08-02 (Date of Last Commit)

* Updated the SnapATAC2 docker to include v2.7.0; the pipeline will now produce a library-level summary metric CSV for the BAM.

# 5.5.0
2024-08-06 (Date of Last Commit)

* Updated the warp-tools docker to calculate mitochondrial reads from unique reads in cell and gene metrics; these metrics are in the cell and gene metrics CSV as well as h5ad

# 5.4.1
2024-08-02 (Date of Last Commit)

* The ubuntu_16_0_4 docker image version was pinned instead of using the latest tag; this does not affect the outputs of the pipeline

# 5.4.0
2024-07-25 (Date of Last Commit)

* Updated the warp-tools docker image to add TSO metrics to the output h5ad and metric CSV files
* Update the library-level metrics to include new TSO metrics and NHashID descriptor

# 5.3.1
2024-07-18 (Date of Last Commit)

* The atac.wdl was refactored into its own directory under the pipelines/skylab directory; this change does not impact the Multiome outputs

# 5.3.0
2024-07-11 (Date of Last Commit)

* Updated the Multiome.wdl to run on Azure; cloud_provider is a new, required input

# 5.2.0
2024-07-09 (Date of Last Commit)

* Added new optional input parameter of nhash_id, an optional identifier for a library aliquot that is echoed in the ATAC fragment h5ad, the gene expression h5ad (in the data.uns), and the gene expression library metrics CSV output; default is set to null
* Added test statements again for GH action (to release from develop). Will probably revert

# 5.1.0
2024-06-28 (Date of Last Commit)

* Updated the STARsolo parameters for estimating cells to Emptydrops_CR
* Added an optional input for expected cells which is used for metric calculation

# 5.0.0
2024-05-20 (Date of Last Commit)

* Updated SnapATAC2 docker to SnapATAC2 v2.6.3; this impacts the workflow output metrics

# 4.0.2
2024-05-14 (Date of Last Commit)

* Updated the Paired-tag Demultiplex task so that some intermediate input names have been renamed; this change does not impact the Multiome workflow

# 4.0.1
2024-05-10 (Date of Last Commit)

* Updated the Paired-tag Demultiplex task; this change does not impact the Multiome workflow

# 4.0.0
2024-04-24 (Date of Last Commit)

* Updated the input parameters for STARsolo in STARsoloFastq task. These include the parameters: soloCBmatchWLtype, soloUMIdedup and soloUMIfiltering
* Added "Uniform" as the default string for STARsolo multimapping parameters

# 3.4.2
2024-04-03 (Date of Last Commit)
* Modified adaptor trimming in Paired-tag WDL; this does not impact Multiome

# 3.4.1
2024-03-26 (Date of Last Commit)

* Updated the median umi per cell metric for STARsolo library-level metrics

# 3.4.0 
2024-03-15 (Date of Last Commit)

* Added cell metrics to the library-level metrics

* Updated the docker for the MergeStarOutput task to include STARsolo v2.7.11a and custom scripts to create a uniform matrix file and scripts to collect library-level metrics from STARsolo output

* Modified the MergeStarOutput to call a custom script for creating a uniform matrix file (mtx) from individual shard mtx files and to create a filtered matrix from the uniform matrix with STARsolo

# 3.3.0
2024-02-28 (Date of Last Commit)

* Added the gene expression library-level metrics CSV as output of the Multiome pipeline; this is produced by the Optimus subworkflow

# 3.2.1
2024-02-29 (Date of Last Commit)

* Moved the disk and mem for the Multiome Join Barcodes task into the task inputs section


# 3.2.0
2024-02-22 (Date of Last Commit)

* Updated StarAlign.MergeStarOutput to add a shard number to the metrics files
* Removed ref_genome_fasta input from Multiome WDL and JSON

# 3.1.3
2024-02-07 (Date of Last Commit)

* Updated the Metrics tasks to exclude mitochondrial genes from reads_mapped_uniquely, reads_mapped_multiple and reads_mapped_exonic, reads_mapped_exonic_as and reads_mapped_intergenic

# 3.1.2
2024-02-01 (Date of Last Commit)

* Add new paired-tag task to parse sample barcodes from cell barcodes when preindexing is set to true; this does not affect the Multiome pipeline


# 3.1.1 
2024-01-30 (Date of Last Commit)

* Added task GetNumSplits before FastqProcess ATAC task to determine the number of splits based on the bwa-mem2 machine specs
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of splits equal the number of ranks
* Added an error message to the BWAPairedEndAlignment ATAC task to ensure that the number of R1s equals to the number of R3s

# 3.1.0
2024-01-24 (Date of Last Commit)
* Promote aligner_metrics from Optimus task level outputs to Multiome pipeline level outputs

# 3.0.5 
2024-01-18 (Date of Last Commit)

* Increased memory for MergeStarOutputs in StarAlign.wdl, RunEmptyDrops in RunEmptyDrops.wdl, OptimusH5ad in H5adUtils.wdl and GeneMetrics in Metrics.wdl
* Added the --soloMultiMappers flag as an optional input to the StarSoloFastq task in the StarAlign.wdl
* Added a check of read2 length to the paired-tag pipeline; this does not affect the Multiome workflow

# 3.0.4
2024-01-05 (Date of Last Commit)

* Added new functionality to the ATAC workflow for paired-tag data, including the option for SnapATAC to pull cell barcodes from the BB tag of the BAM
* Modified the STARsoloFastq task in the StarAlign.wdl so STARsolo can run different types of alignments in a single STARsolo command depending on the counting_mode

# 3.0.3
2023-12-20 (Date of Last Commit)

* Added updated docker to BWAPairedEndAlignment ATAC task to use updated code for distributed bwa-mem2 from Intel
* Removed MergedBAM ATAC and moved BWAPairedEndAlignment ATAC outside of the for loop
* Changed CPU platform to Ice Lake for BWAPairedEndAlignment ATAC task
* Added input parameter input_output_parameter to the Multiome ATAC wdl

# 3.0.2
2023-12-20 (Date of Last Commit)

* JoinMultiomeBarcodes now has dynamic memory and disk allocation

# 3.0.1
2023-12-12 (Date of Last Commit)

* ValidateVcfs now has optional memory parameter; this does not affect this pipeline
* Downgraded Cell Bender from v0.3.1 to v0.3.0

# 3.0.0
2023-11-22 (Date of Last Commit)

* Added gene expression barcodes to the Multiome ATAC fragment file
* Updated the JoinBarcodes task to bgzip and tabix the final ATAC fragment file
* Added the tabix index file as an output to Multiome

# 2.3.3
2023-11-21 (Date of Last Commit)

* Added the latest warp-tools docker to tasks in the Metrics, FastqProcessing and H5adUtils wdls; this incorporates new input parameter for number of output fastq files to fastqprocess

# 2.3.2
2023-11-20 (Date of Last Commit)

* Added an optional task to the Multiome.wdl that will run CellBender on the Optimus output h5ad file

# 2.3.1
2023-11-20 (Date of Last Commit)

* Added the latest warp-tools docker to the Metrics task; this allows use of REFSEQ references

# 2.3.0
2023-11-03 (Date of Last Commit)

* Updated the Metrics task so that Cell Metrics and Gene Metrics now calculate intronic, intronic_as, exonic, exonic_as, and intergenic metrics from unique reads only using the NH:i:1 tag in the BAM

# 2.2.2
2023-10-20 (Date of Last Commit)

* Removed Dropna from JoinBarcodes subtask of the H5adUtils task, which was causing the JoinBarcodes to fail for some gene expression matrices

* Updated path to Multiome whitelists to reflect location in public storage.

# 2.2.0
2023-10-05 (Date of Last Commit)
* Added a JoinMultiomeBarcodes task to the H5adUtils that adds a column in the ATAC and Optimus output h5ad linking gene expression and ATAC barcodes

# 2.1.0
2023-09-21 (Date of Last Commit)
* Added dynamic barcode orientation selection to the ATAC workflow FastqProcess task

# 2.0.0
2023-09-05 (Date of Last Commit)

* Updated Optimus pipeline to include STARsolo v2.7.11a
* Added sF tag to STARsolo aligner parameters
* Updated TagSort tool for Optimus Metrics task to calculate metrics based on the sF tag
* Modified H5adUtils task to include new metrics in the final Optimus h5ad
* Removed the Dropseq metrics task
* Updated the ATAC wdl to optimize the MakeFragmentFile task for cpu, memory, disk, and cpu platform
 
# 1.0.1 
2023-07-23 (Date of Last Commit)

* Added STARsolo v2.7.10b metric outputs as an optional pipeline output and an output of the STARalign and MergeSTAR tasks

* Updated the CountAlignments task in the FeatureCounts.wdl to use a new docker image. This change does not affect the Multiome pipeline

# 1.0.0
2023-06-22 (Date of Last Commit)

* Initial release of the multiome pipeline

