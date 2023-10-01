version 1.0

workflow MyGetChromosomes {
  input {
    Array[String] unique_chromosomes 
    Array[String] chromosomes
    Array[String] ordered_vcf_files
  }

  scatter (chrom in unique_chromosomes) {
    scatter (p in zip(ordered_vcf_files, chromosomes)) {
      if (p.right == chrom) {
        String vcf = p.left
      } 
    }
    Pair[String,Array[String]] group = (chrom, select_all(vcf))
  }

  output {
    Array[Pair[String,Array[String]]] chrom_vcf_map = group
  }
}
