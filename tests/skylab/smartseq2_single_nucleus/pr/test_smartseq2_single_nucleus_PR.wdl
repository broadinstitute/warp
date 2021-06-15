version 1.0

import "../../../../pipelines/skylab/smartseq2_single_nucleus/SmartSeq2SingleNucleus.wdl" as target_wdl
import "../../../../tests/skylab/smartseq2_single_nucleus/pr/ValidateSmartSeq2SingleNucleus.wdl" as checker_wdl
import "../../../../verification/VerifyTasks.wdl" as verify_tasks

# this task will be run by the jenkins script that gets executed on our PRs.
workflow TestSmartSeq2SingleNucleusPR {
  input {

    #checksums
    String truth_exon_intron_counts_hash

    File truth_bam
    File truth_loom


    # snSS2 inputs
    File genome_ref_fasta
    File star_reference
    File annotations_gtf
    String input_id
    String output_name
    File adapter_list
    File fastq1
    File fastq2
  }

  call target_wdl.SmartSeq2SingleNucleus as target_workflow {
    input:
      genome_ref_fasta = genome_ref_fasta,
      input_id = input_id,
      output_name = output_name,
      fastq1 = fastq1,
      fastq2 = fastq2,
      adapter_list = adapter_list,
      annotations_gtf = annotations_gtf,
      star_reference = star_reference

  }

   call checker_wdl.ValidateSnSmartSeq2 as checker_workflow {
     input:
      exon_intron_counts_hash = target_workflow.exon_intron_counts,
      truth_exon_intron_counts_hash = truth_exon_intron_counts_hash,
      loom_output = target_workflow.loom_output_file,
      truth_loom = truth_loom
   }

   call verify_tasks.CompareBams as CompareBams {
     input:
       test_bam = target_workflow.aligned_bam,
       truth_bam = truth_bam
     }

}
