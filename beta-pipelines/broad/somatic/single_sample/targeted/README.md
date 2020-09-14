# TargetedSomaticSingleSample Pipeline
This WDL pipeline implements data pre-processing according to the GATK Best Practices 
for somatic human targeted sequencing data.

## Requirements/Expectations :
- Human targeted sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
  - filenames all have the same suffix (we use ".unmapped.bam")
  - files must pass validation by ValidateSamFile
  - reads are provided in query-sorted order
  - all reads must have an RG tag
- Reference genome must be Hg38 with ALT contigs

Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
For program versions, see docker containers.

## Default Input Values :
- By default ApplyBQSR will NOT bin the base qualities for the targeted somatic pipeline.
  - To override this, add '"TargetedSomaticSingleSample.bin_base_qualities": true' to the inputs.json file

- By default reads will NOT be hardclipped to remove adapter sequence. Hard Clipping is recommended for TWIST somatic exome data.
  - To override this, add '"TargetedSomaticSingleSample.hard_clip_reads": true' to the inputs.json file
