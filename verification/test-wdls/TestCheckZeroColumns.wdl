version 1.0

import "../VerifyOptimus.wdl" as VerifyOptimus

workflow TestCheckZeroColumns {
  input {
    File gene_metrics
    File cell_metrics
    File library_metrics
  }

  call VerifyOptimus.CheckForZeroColumns as CheckZero {
    input:
      gene_metrics = gene_metrics,
      cell_metrics = cell_metrics,
      library_metrics = library_metrics
  }

  output {
    String result = CheckZero.result
  }
}
