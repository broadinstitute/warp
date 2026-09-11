# MTtoVCF
The MTtoVCF pipelines take a specific hail matrix table, and does the following.
The specific matrix table is referred to as the superset matrix table, consisting of 10750 multiomic samples.
Ths MTtoVCF.wdl and ComputeGeneticPCs.wdl should be generic to other matrix tables.
1. FilterMT.wdl
  - Filters the matrix down to specified sample list
  - Filter the matrix to positions where AC > some given value
  - Filter on some basic quality criteria
  - Outputs a checkpoint matrix table to a specified location
  - This breaks up the processing of the matrix table into two steps to allow us to take advantage of checkpointing in case of failures, and also allows for a faster read/write in the following step.
2. MTtoVCF
  - Take a matrix table and outputs a VCF representation to a specified location.
3. ComputeGeneticPCs.wdl
  - Takes a VCF and uses plink2? to compute genetic pcs. Need to verify if plink2 was the correct tool for this.
