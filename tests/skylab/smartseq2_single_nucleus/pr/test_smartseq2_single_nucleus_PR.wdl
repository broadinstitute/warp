version 1.0

import https://raw.githubusercontent.com/broadinstitute/warp/snSS2_first_wdls/pipelines/skylab/smartseq2_single_nucleus/SmartSeq2SingleNucleus.wdl as target_wdl
import "../../../../tests/skylab/smartseq2_single_nucleus/pr/ValidateSmartSeq2SingleNucleus.wdl" as checker_wdl


# this task will be run by the jenkins script that gets executed on our PRs.
workflow TestSmartSeq2SingleNucleusPR {
  input {
    # Validation input
    #String loom_output
    File counts
    String expected_counts_hash
    File? target_metrics
    String expected_metrics_hash

    # snSS2 inputs
    File genome_ref_fasta
    File tar_star_reference
    File annotations_gtf
    String stranded
    String input_id
    String output_name
    File adapter_list
    File fastq1
    File? fastq2
    Boolean paired_end
  }

  call target_wdl.SmartSeq2SingleNucleus as target_workflow {
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

  call checker_wdl.ValidateSmartSeq2Plate as checker_workflow {
    input:
      loom_output = target_workflow.loom_output,
      truth_loom = truth_loom
  }

  call checker_wdl.ValidateSmartSeq2Plate as checker_workflow {
      input:
        loom_output = target_workflow.loom_output,
        truth_loom = truth_loom
    }
}