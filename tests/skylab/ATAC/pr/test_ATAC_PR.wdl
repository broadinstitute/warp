version 1.0

import "../../../../beta-pipelines/skylab/ATAC/ATAC.wdl" as target
import "../../../../tests/skylab/ATAC/pr/ValidateATAC.wdl" as checker

# this workflow will be run by the jenkins script that gets executed by PRs.
workflow TestAtacPR {
  input {
      # output hashes
      String expected_snap_hash
      String expected_snapqc_hash
      String expected_bam_hash

      # scATAC inputs
      File input_fastq1
      File input_fastq2
      String bin_size_list

      # Reference related variables
      File input_reference
      String genome_size_file
      String genome_name

      # outputbam not used at the minute
      String output_bam

  }

  call target.ATAC as target {
    input:
       fastq_gzipped_input_read1 = input_fastq1,
       fastq_gzipped_input_read2 = input_fastq2,
       bin_size_list = bin_size_list,
       tar_bwa_reference = input_reference,
       min_map_quality = 30,
       quality_cutoff = 0,
       adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
       adapter_seq_read2 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
       min_length = 10,
       output_base_name = "scATAC",
       tar_bwa_reference = input_reference,
       genome_size_file = genome_size_file,
       genome_name = genome_name,
       max_fragment_length = 2000
  }

  call checker.ValidateATAC as checker {
    input:
        snap = target.snap_output,
        snapqc = target.snap_qc_output,
        bam = target.bam_filtered_and_sorted_compliant_output,
        ## Todo: add mitochondiral bam
        expected_snap_hash = expected_snap_hash,
        expected_snapqc_hash = expected_snapqc_hash,
        expected_bam_hash = expected_bam_hash
  }

}
