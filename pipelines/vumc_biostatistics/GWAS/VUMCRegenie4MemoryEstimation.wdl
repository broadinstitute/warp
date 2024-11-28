version 1.0

import "./GWASUtils.wdl" as GWASUtils

workflow VUMCRegenie4MemoryEstimation {
  input {
    Int num_sample
    Int num_variant
    Int num_phenotype
    Int num_chromosome
    Int num_covariate
    Int num_ridge_l0 = 5
    Int block_size = 1000
  }

  call GWASUtils.Regenie4MemoryEstimation {
    input:
      num_sample = num_sample,
      num_variant = num_variant,
      num_phenotype = num_phenotype,
      num_chromosome = num_chromosome,
      num_covariate = num_covariate,
      num_ridge_l0 = num_ridge_l0,
      block_size = block_size
  }

  output {
    Int step1_memory_gb = Regenie4MemoryEstimation.step1_memory_gb
  }
}
