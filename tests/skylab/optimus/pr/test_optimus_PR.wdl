version 1.0

import "../../../../pipelines/skylab/optimus/Optimus.wdl" as target
import "../../../../tests/skylab/optimus/pr/ValidateOptimus.wdl" as checker


# this workflow will be run by the jenkins script that gets executed by PRs.
workflow TestOptimusPR {
  input {
    # output hashes
    String expected_bam_hash
    String expected_gene_metric_hash
    String expected_cell_metric_hash
    String expected_loom_file_checksum
    File reference_matrix

    # Optimus inputs
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq

    File whitelist  # 10x genomics cell barcode whitelist for 10x V2
    File tar_star_reference  # star reference
    File annotations_gtf  # gtf containing annotations for gene tagging
    File ref_genome_fasta  # genome fasta file
    String input_id  # name of sample matching this file, inserted into read group header
    String chemistry # chemistry identifier

    Boolean force_no_check
  }
  

  call target.Optimus as target {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      i1_fastq = i1_fastq,
      whitelist = whitelist,
      tar_star_reference = tar_star_reference,
      annotations_gtf = annotations_gtf,
      ref_genome_fasta = ref_genome_fasta,
      input_id = input_id,
      chemistry = chemistry,
      force_no_check = force_no_check,
      emptydrops_lower = 1
  }

  call checker.ValidateOptimus as checker {
    input:
      bam = target.bam,
      matrix = target.matrix,
      matrix_row_index = target.matrix_row_index,
      matrix_col_index = target.matrix_col_index,
      gene_metrics = target.gene_metrics,
      cell_metrics = target.cell_metrics,
      loom_file = target.loom_output_file,
      
      reference_matrix = reference_matrix,
      expected_bam_hash = expected_bam_hash,
      expected_cell_metric_hash = expected_cell_metric_hash,
      expected_gene_metric_hash = expected_gene_metric_hash,
      expected_loom_file_checksum = expected_loom_file_checksum
  }

}
