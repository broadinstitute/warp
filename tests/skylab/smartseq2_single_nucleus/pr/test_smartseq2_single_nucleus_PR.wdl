version 1.0

import "../../../../pipelines/skylab/smartseq2_single_nucleus_multisample/MultiSampleSmartSeq2SingleNucleus.wdl" as target_wdl
import "../../../../tests/skylab/smartseq2_single_nucleus/pr/ValidateSmartSeq2SingleNucleus.wdl" as checker_wdl
import "../../../../verification/VerifyTasks.wdl" as verify_tasks

# this task will be run by the jenkins script that gets executed on our PRs.
workflow TestSmartSeq2SingleNucleusPR {
  input {

    #checksums
    Array[String] truth_exon_intron_counts_hash
    Array[File] truth_bam
    File truth_loom


    # snSS2 inputs
    File genome_ref_fasta
    File tar_star_reference
    File annotations_gtf
    File adapter_list
    String batch_id
    Array[String] input_ids
    Array[File] fastq1_input_files
    Array[File] fastq2_input_files
  }

  call target_wdl.MultiSampleSmartSeq2SingleNucleus as target_workflow {
    input:
      genome_ref_fasta = genome_ref_fasta,
      input_ids = input_ids,
      batch_id = batch_id,
      fastq1_input_files = fastq1_input_files,
      fastq2_input_files = fastq2_input_files,
      adapter_list = adapter_list,
      annotations_gtf = annotations_gtf,
      tar_star_reference = tar_star_reference
  }
  scatter(idx in range(length(input_ids))) {
    call checker_wdl.ValidateSnSmartSeq2 as checker_workflow {
      input:
        exon_intron_counts_hash = target_workflow.exon_intron_count_files[idx],
        truth_exon_intron_counts_hash = truth_exon_intron_counts_hash[idx],
        loom_output = target_workflow.loom_output,
        truth_loom = truth_loom
    }
    call verify_tasks.CompareBams as CompareBams {
      input:
        test_bam = target_workflow.bam_files[idx],
        truth_bam = truth_bam[idx],
        lenient_header = true
    }
    }

}
