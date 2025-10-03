# snM3C Test Data
The following lists the data sets used for workflow validation, including workflow engineering (plumbing) tests and scientific validation tests.

* WARP tests for snM3C use data downsampled to two cells (M16 barcode and G13 barcode) using the Utility WDL [sample_fastq](../../../../tasks/skylab/sample_fastq.14.wdl). 

  * The whitelist used for the downsampling workflow is "gs://fc-f62be246-9bf6-42a8-9123-7bd1c52894a7/downsampled/M16_G13_whitelist.txt", which is simply a list of the M16 and G13 barcodes.
  * The read structure for downsample is "8C148M."

* To ensure tests run quickly, the barcode list has been subsetted to the barcodes for these two cells.

* Both original data sets (FASTQs) were provided courtesy of Dr. Hanqing Liu.

* Original raw files are privately available in [this Terra worksapce](). 

## Plumbing data
The plumbing data is downsampled from a human multiplexed (384-well plate) Miseq sample called 230117-AD3C-hs-snm3C_seq-NovaSeq-pe-150-PB-AD3C_BA17_2027_P1-1-B11. The FASTQs represent four lanes of sequencing. 

## Scientific data
The scientific data is downsampled from a human multiplexed (384-well plate) Novaseq sample called 230419-iN-hs-snm3C_seq-NovaSeq-pe-150-BW-iN230412_Entorhinal_Cortex_Adult_2115_R1_1-1-I3. The FASTQs represent four lanes of sequencing. 
