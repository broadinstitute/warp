version 1.0

import "../../../../pipelines/skylab/smartseq2_single_sample/SmartSeq2SingleSample.wdl" as target_wdl
import "../../../../tests/skylab/smartseq2_single_sample/pr/ValidateSmartSeq2SingleCell.wdl" as checker_wdl

# this task will be run by the jenkins script that gets executed on our PRs.
workflow TestSmartSeq2SingleCellPR {
  input { 
    # expected hashes of target_workflow outputs
    String expected_counts_hash
    String expected_metrics_hash = ""

    # SS2 inputs
    File genome_ref_fasta
    File rrna_intervals
    File gene_ref_flat
    File hisat2_ref_index
    File hisat2_ref_trans_index
    File rsem_ref_index
    String hisat2_ref_name
    String hisat2_ref_trans_name
    String stranded
    String input_id
    String output_name
    File fastq1
    File fastq2
  }

  call target_wdl.SmartSeq2SingleCell as target_workflow {
    input:
      genome_ref_fasta = genome_ref_fasta,
      rrna_intervals = rrna_intervals,
      gene_ref_flat = gene_ref_flat,
      hisat2_ref_index = hisat2_ref_index,
      hisat2_ref_trans_index = hisat2_ref_trans_index,
      rsem_ref_index = rsem_ref_index,
      hisat2_ref_name = hisat2_ref_name,
      hisat2_ref_trans_name = hisat2_ref_trans_name,
      stranded = stranded,
      input_id = input_id,
      output_name = output_name,
      fastq1 = fastq1,
      fastq2 = fastq2,
      paired_end = true
  }

  call checker_wdl.ValidateSmartSeq2SingleCell as checker_workflow {
    input:
     counts = target_workflow.rsem_gene_results,
     expected_counts_hash = expected_counts_hash,
     target_metrics = target_workflow.insert_size_metrics,
     expected_metrics_hash = expected_metrics_hash
  }

}
