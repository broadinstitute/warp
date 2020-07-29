version 1.0

import "MultiSampleSmartSeq2.wdl" as target_wdl
import "ValidateMultiSampleSmartSeq2.wdl" as checker_wdl

workflow TestMultiSampleSmartSeq2 {
  input {
      # SS2 inputs
      File genome_ref_fasta
      File rrna_intervals
      File gene_ref_flat
      String hisat2_ref_name
      String hisat2_ref_trans_name
      File hisat2_ref_index
      File hisat2_ref_trans_index
      File rsem_ref_index
      String stranded
      Boolean paired_end
      File truth_loom

      # Plate information and input files
      String file_prefix
      Array[String] input_file_names
      String batch_id
  }

  call target_wdl.MultiSampleSmartSeq2 as target_workflow {
    input:
      genome_ref_fasta = genome_ref_fasta,
      rrna_intervals = rrna_intervals,
      gene_ref_flat = gene_ref_flat,
      hisat2_ref_name = hisat2_ref_name,
      hisat2_ref_trans_name = hisat2_ref_trans_name,
      hisat2_ref_index = hisat2_ref_index,
      hisat2_ref_trans_index = hisat2_ref_trans_index,
      rsem_ref_index = rsem_ref_index,
      stranded = stranded,
      file_prefix = file_prefix,
      input_file_names = input_file_names,
      batch_id = batch_id,
      paired_end = paired_end
  }

  call checker_wdl.ValidateSmartSeq2Plate as checker_workflow {
    input:
      loom_output = target_workflow.loom_output,
      truth_loom = truth_loom
  }
}
