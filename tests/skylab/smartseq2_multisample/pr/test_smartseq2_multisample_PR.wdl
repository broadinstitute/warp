version 1.0

import "../../../../pipelines/skylab/smartseq2_multisample/MultiSampleSmartSeq2.wdl" as target_wdl
import "../../../../tests/skylab/smartseq2_multisample/pr/ValidateMultiSampleSmartSeq2.wdl" as checker_wdl

workflow TestMultiSampleSmartSeq2 {
  input {
        # Gene Annotation
        File genome_ref_fasta
        File rrna_intervals
        File gene_ref_flat

        # Reference index information
        File hisat2_ref_name
        File hisat2_ref_trans_name
        File hisat2_ref_index
        File hisat2_ref_trans_index
        File rsem_ref_index

        # Sample information
        String stranded
        Array[String] input_ids
        Array[String]? input_names
        Array[String] fastq1_input_files
        Array[String] fastq2_input_files = []
        String batch_id
        String? batch_name
        String? input_name_metadata_field
        String? input_id_metadata_field
        Boolean paired_end

        # Validation input
        File truth_loom
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
      input_ids = input_ids,
      input_names = input_names,
      fastq1_input_files = fastq1_input_files,
      fastq2_input_files = fastq2_input_files,
      batch_id = batch_id,
      batch_name = batch_name,
      input_name_metadata_field = input_name_metadata_field,
      input_id_metadata_field = input_id_metadata_field,
      paired_end = paired_end
  }

  call checker_wdl.ValidateSmartSeq2Plate as checker_workflow {
    input:
      loom_output = target_workflow.loom_output,
      truth_loom = truth_loom
  }
}
