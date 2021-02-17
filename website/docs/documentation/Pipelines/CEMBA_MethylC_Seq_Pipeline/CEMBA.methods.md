# CEMBA_v1.1.0 Publication Methods

Below we provide a sample methods section for a publication. For the complete pipeline documentation, see the [CEMBA README](./README.md).

## Methods

Data processing was performed with the CEMBA v1.1.0 Pipeline. Sequencing reads were first trimmed to remove adaptors using Cutadapt 1.18 with the following parameters in paired-end mode: `-f fastq -quality-cutoff 20 -minimum-length 62 -a AGATCGGAAGAGCACACGTCTGAAC -A AGATCGGAAGAGCGTCGTGTAGGGA`.

After trimming the adapters, an unaligned BAM (uBAM) for the trimmed R1 FASTQ was created using Picard v2.18.23.

Cell barcodes were then extracted from the trimmed R1 FASTQ and tagged to the R1 uBAM with Single Cell Tools (sctools) v0.3.4a using a barcode whitelist as well as configurable barcode start positions and lengths.

Next, for multiplexed samples, the random primer index sequence and Adaptase C/T tail were further removed from the adaptor-trimmed R1 and R2 FASTQs using Cutadapt with the following parameters: `-f fastq -quality-cutoff 16 -quality-cutoff -16 -minimum-length 30`.

The trimmed R1 and R2 reads were then aligned to mouse (mm10) or human (hg19) genomes separately as single-end reads using Bismark v0.21.0 with the parameters `--bowtie2 --icpc --X 2000` (paired-end mode) and `--pbat` (activated for mapping R1 reads).

After alignment, the output R1 and R2 BAMs were sorted in coordinate order and duplicates removed using the Picard MarkDuplicates REMOVE_DUPLICATE option. Samtools 1.9 was used to further filter BAMs with a minimum map quality of 30 using the parameter `-bhq 30`.

Methylation reports were produced for the filtered BAMs using Bismark. The barcodes from the R1 uBAM were then attached to the aligned, filtered R1 BAM with Picard. The R1 and R2 BAMs were merged with Samtools. Readnames were added to the merged BAM and a methylated VCF created using MethylationTypeCaller in GATK 4.1.2.0. The VCF was then converted to an additional ALLC file using a custom python script. 

Samtools was then used to calculate coverage depth for sites with coverage greater than 1 and to create BAM index files. The final outputs included the barcoded aligned BAM, BAM index, a VCF with locus-specific methylation information, VCF index, ALLC file, and methylation reports.

An example of the pipeline and its outputs is available on [Terra](https://app.terra.bio/#workspaces/brain-initiative-bcdc/Methyl-c-seq_Pipeline). Examples of genomic reference files and other inputs can be found in the pipeline’s [example JSON](https://github.com/broadinstitute/warp/blob/develop/pipelines/cemba/cemba_methylcseq/example_inputs/CEMBA.inputs.json).
