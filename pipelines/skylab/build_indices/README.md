# Building References for Single-Cell Pipelines

The macaque and mouse reference files for the Optimus pipeline (used by BICAN and HCA) can be prepared by using 
the `BuildIndices.wdl` workflow. The WDL is self contained without any dependent WDL files. The references can be built
 by providing an input JSON file to the above WDL.

Here are two example input files for the above WDL:

`Macaque.json`:
```json
{
 "BuildIndices.annotations_gtf":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Macaque/genomic.gtf",
 "BuildIndices.biotypes":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Biotypes.tsv",
 "BuildIndices.genome_fa":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Macaque/GCF_003339765.1_Mmul_10_genomic.fna",
 "BuildIndices.gtf_annotation_version":"103",
 "BuildIndices.genome_source":"NCBI",
 "BuildIndices.genome_build":"GCF_003339765.1",
 "BuildIndices.organism":"Macaque"
}

```

`Mouse.json`:
```json
{
 "BuildIndices.annotations_gtf":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Mouse/gencode.vM32.primary_assembly.annotation.gtf",
 "BuildIndices.biotypes":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Biotypes.tsv",
 "BuildIndices.genome_fa":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Mouse/GRCm39.primary_assembly.genome.fa",
 "BuildIndices.gtf_annotation_version":"M32",
 "BuildIndices.genome_source":"GENCODE",
 "BuildIndices.genome_build":"GRCm39",
 "BuildIndices.organism":"Mouse"
}


```

The workflow has one task:

- `BuildStarSingleNucleus` - Builds the reference for STAR aligner and creates a modified GTF for the Single Nucleus Smart-Seq2 Pipeline.
