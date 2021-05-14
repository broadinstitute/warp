version 1.0

import "https://raw.githubusercontent.com/broadinstitute/warp/np_snss2_pr_test/pipelines/skylab/smartseq2_single_nucleus/SmartSeq2SingleNucleus.wdl" as target_wdl
import "https://raw.githubusercontent.com/broadinstitute/warp/np_snss2_pr_test/tests/skylab/smartseq2_single_nucleus/pr/ValidateSmartSeq2SingleNucleus.wdl" as checker_wdl
import "https://raw.githubusercontent.com/broadinstitute/warp/develop/verification/VerifyTasks.wdl" as verify_tasks
import "https://raw.githubusercontent.com/broadinstitute/warp/develop/tests/skylab/smartseq2_multisample/pr/ValidateMultiSampleSmartSeq2.wdl" as validate

# this task will be run by the jenkins script that gets executed on our PRs.
workflow TestSmartSeq2SingleNucleusPR {
  input {

    #checksums
    String introns_counts_hash
    String expected_introns_counts_hash
    String exons_counts_hash
    String expected_exons_counts_hash

    File truth_bam
    File truth_loom


    # snSS2 inputs
    File genome_ref_fasta
    File tar_star_reference
    File annotations_gtf
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
      input_id = input_id,
      output_name = output_name,
      fastq1 = fastq1,
      fastq2 = fastq2,
      paired_end = true,
      adapter_list = adapter_list,
      annotations_gtf = annotations_gtf,
      tar_star_reference = tar_star_reference

  }

  call validate.ValidateSmartSeq2Plate as compareLooms {
    input:
      loom_output = target_workflow.loom_output_files,
      truth_loom = truth_loom
  }

   call checker_wdl.CompareCounts as checker_workflow {
     input:
      introns_counts_hash = introns_counts_hash,
      expected_introns_counts_hash = expected_introns_counts_hash,
      exons_counts_hash = exons_counts_hash,
      expected_exons_counts_hash = expected_exons_counts_hash,
      loom_output = target_workflow.loom_output_files,
      truth_loom = truth_loom
      # target_metrics = target_workflow.insert_size_metrics,
      # expected_metrics_hash = expected_metrics_hash
   }

   call verify_tasks.CompareBams as CompareBams {
     input:
       test_bam = target_workflow.aligned_bam,
       truth_bam = truth_bam
     }

}