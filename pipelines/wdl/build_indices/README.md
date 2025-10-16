# Building References for Single-Cell Pipelines

The macaque and mouse reference files for the Optimus pipeline (used by BICAN and HCA) can be prepared by using 
the `BuildIndices.wdl` workflow. The WDL is self contained without any dependent WDL files. The references can be built
 by providing an input JSON file to the above WDL.

Here are two example input files for the above WDL:

`Macaque.json`:
```json
{
 "BuildIndices.annotations_gtf":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Macaque/genomic.gtf",
 "BuildIndices.biotypes":"gs://fc-df68cb43-8c48-401b-9ef1-7cbb3acc788d/Biotypes.tsv",
 "BuildIndices.genome_fa":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Macaque/GCF_003339765.1_Mmul_10_genomic.fna",
 "BuildIndices.gtf_version":"GCF_003339765.1",
 "BuildIndices.organism":"Macaque",
 "BuildIndices.organism_prefix":"Macaque"
}
```

`Mouse.json`:
```json
{
 "BuildIndices.BuildStarSingleNucleus.organism":"Mouse",
 "BuildIndices.BuildStarSingleNucleus.organism_prefix":"Mouse",
 "BuildIndices.annotations_gtf":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Mouse/gencode.vM31.primary_assembly.annotation.gtf",
 "BuildIndices.biotypes":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Biotypes.tsv",
 "BuildIndices.genome_fa":"gs://fc-c40ec8a8-d60f-42f7-be36-3986b475190a/Mouse/GRCm39.primary_assembly.genome.fa",
 "BuildIndices.gtf_version":"M31"
}

```

The workflow has one task:

- `BuildStarSingleNucleus` - Builds the reference for STAR aligner and creates a modified GTF for the Single Nucleus Smart-Seq2 Pipeline.
