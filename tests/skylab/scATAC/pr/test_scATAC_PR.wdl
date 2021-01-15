version 1.0

import "../../../../pipelines/skylab/scATAC/scATAC.wdl" as target
import "../../../../tests/skylab/scATAC/pr/ValidateSCATAC.wdl" as checker

# this workflow will be run by the jenkins script that gets executed by PRs.
workflow TestScAtacPR {
  input {
      # output hashes
      String expected_snap_hash
      String expected_snapqc_hash
      String expected_bam_hash

      # scATAC inputs
      File input_fastq1
      File input_fastq2
      File input_reference
      String output_bam
      String genome_name
      String bin_size_list
      String input_id
  }

  call target.scATAC as target {
    input:
        input_fastq1 = input_fastq1,
        input_fastq2 = input_fastq2,
        input_reference = input_reference,
        output_bam = output_bam,
        genome_name = genome_name,
        bin_size_list = bin_size_list,
        input_id = input_id
  }

  call checker.ValidateSCATAC as checker {
    input:
        snap = target.output_snap,
        snapqc = target.output_snap_qc,
        bam = target.output_aligned_bam,

        expected_snap_hash = expected_snap_hash,
        expected_snapqc_hash = expected_snapqc_hash,
        expected_bam_hash = expected_bam_hash
  }

}
