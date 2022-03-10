# Building References for HCA Pipelines

The human and mouse reference files for the HCA pipelines (HCA SS2 and HCA Optimus pipeline) can be prepared by using 
the `BuildIndices.wdl` workflow. The WDL is self contained without any dependent WDL files. The references can be built
 by providing an input JSON file to the above WDL.

Here are two example input files for the above WDL:

`mouse_inputs.json`:
```json
{
  "BuildIndices.gtf_version": "M21",
  "BuildIndices.organism": "mouse",
  "BuildIndices.organism_prefix": "m",
  "BuildIndices.genome_short_string": "mm10",
  "BuildIndices.dbsnp_version": "150"
}
```

`human_inputs.json`:
```json
{
  "BuildIndices.gtf_version": "27",
  "BuildIndices.organism": "human",
  "BuildIndices.organism_prefix": "h",
  "BuildIndices.genome_short_string": "hg38",
  "BuildIndices.dbsnp_version": "150"
}
```

The workflow has several tasks:

- `GetReference` - This task fetches the reference files: the annotation file and the genome file.

- `BuildPicardRefFlat` - This step creates a flat file from the annotation file for Picard. (HCA SS2)

- `BuilidHisat2SnpHaplotype` - Builds the HISAT2 reference with the SNP information. (HCA SS2)

- `BuildIntervalList` - This step creates a text files with intervals of the chromosomes to be used with Picard. (HCA SS2)

- `BuildHisat2` - Builds the HISAT2 references from the genomics reference. (HCA SS2)

- `BuildRsem` - This step builds the RSEM references from the annotation and the genome file along with a Bowtie index.
At this point, the Bowtie index is not used. (HCA SS2)

- `BuildStar` - Builds the STAR reference for the Optimus pipeline. (HCA Optimus)

- `BuildHisat2FromRsem` - Builds the reference for HISAT2 from the rsem transcripts sequences. (HCA SS2)

- `BuildStarSingleNucleus` - Builds the reference for STAR aligner and creates a modified GTF for the Single Nucleus Smart-Seq2 Pipeline.
