# WDL Input Overview


This workflow splits a multi-sample VCF file into individual single-sample VCF files. This is achieved by first 
extracting all the samples from the provided multi-sample VCF file. The extracted samples are then broken up into 
"chunks" based on the specified chunk size. A new sample list is generated per chunk. Finally, this new "chunked" 
sample list is used to extract the corresponding samples from the original multi-sample VCF file, creating a new 
single-sample VCF files. These files (and corresponding index files, if requested) are then copied to the specified 
output location. 

## Inputs:

| Input Name           | Description                                                                                                                                                      | Type    | Required | Default                                                                     |
|----------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|----------|-----------------------------------------------------------------------------|
| **multiSampleVcf**   | The GCP path to the multi-sample VCF.                                                                                                                            | File    | Yes      | N/A                                                                         |
| **createIndexFiles** | Whether index files should be generated for each extracted single-sample VCF file.                                                                               | Boolean | Yes      | N/A                                                                         |
| **outputLocation**   | Where the single-sample VCFs (and indices, if requested) should be written to. Note, the user running the workflow must have write permissions to this location. | String  | Yes      | N/A                                                                         |
| **chunkSize**        | The number of samples that should be in each chunk.                                                                                                              | Int     | No       | 1000                                                                        |
| **cpu**              | The number of CPU cores to allocate.                                                                                                                             | Int     | No       | 1                                                                           |
| **memoryMb**         | Memory allocation in megabytes.                                                                                                                                  | Int     | No       | 6000                                                                        |
| **bcftoolsDocker**   | Docker image containing bcftools.                                                                                                                                | String  | No       | us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889 |


## Outputs:
This workflow does not generate any direct outputs. The generated output VCFs (and indices, if requested) will be 
written to the provided output location. The logs, which include any errors or warnings, can be reviewed in the stderr 
file for additional details about the process.