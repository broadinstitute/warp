version 1.0

import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../pipelines/skylab/optimus/Optimus.wdl" as Optimus
import "../../verification/VerifyOptimus.wdl" as VerifyOptimus
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestOptimus {

  input {
    
    # Mode for counting either "sc_rna" or "sn_rna"
    String counting_mode = "sc_rna"

    # Sequencing data inputs
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    String input_id
    String output_bam_basename = input_id
    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field
    # organism reference parameters
    File tar_star_reference
    File annotations_gtf
    File ref_genome_fasta
    File? mt_genes

    # Chemistry options include: 2 or 3
    Int tenx_chemistry_version = 2
    # Whitelist is selected based on the input tenx_chemistry_version

    # Emptydrops lower cutoff
    Int emptydrops_lower = 100

    # Set to true to override input checks and allow pipeline to proceed with invalid input
    Boolean force_no_check = false
    
    # Check that tenx_chemistry_version matches the length of the read 1 fastq;
    # Set to true if you expect that r1_read_length does not match length of UMIs/barcodes for 10x chemistry v2 (26 bp) or v3 (28 bp).
    Boolean ignore_r1_read_length = false

    # Set to true to count reads in stranded mode
    String use_strand_info = "false"
    
# Set to true to count reads aligned to exonic regions in sn_rna mode
    Boolean count_exons = false

    # this pipeline does not set any preemptible varibles and only relies on the task-level preemptible settings
    # you could override the tasklevel preemptible settings by passing it as one of the workflows inputs
    # for example: `"Optimus.StarAlign.preemptible": 3` will let the StarAlign task, which by default disables the
    # usage of preemptible machines, attempt to request for preemptible instance up to 3 times. 

    # Injected from test framework
    String truth_path
    String results_path
    Boolean update_truth
    String vault_token_path
    String google_account_vault_path

  }

  meta {
    allowNestedInputs: true
  }

  call Optimus.Optimus {
    input:
      counting_mode              = counting_mode,
      r1_fastq                   = r1_fastq,
      r2_fastq                   = r2_fastq,
      i1_fastq                   = i1_fastq,
      input_id                   = input_id,
      output_bam_basename        = output_bam_basename,
      input_name                 = input_name,
      input_id_metadata_field    = input_id_metadata_field,
      input_name_metadata_field  = input_name_metadata_field,
      tar_star_reference         = tar_star_reference,
      annotations_gtf            = annotations_gtf,
      ref_genome_fasta           = ref_genome_fasta,
      tenx_chemistry_version     = tenx_chemistry_version,
      emptydrops_lower           = emptydrops_lower,
      force_no_check             = force_no_check,
      use_strand_info            = use_strand_info,
      count_exons                = count_exons,
      ignore_r1_read_length      = ignore_r1_read_length,
  }

  # Collect all of the pipeling output into single Array
  Array[String] pipeline_outputs = select_all([  
                                      Optimus.bam,
                                      Optimus.matrix,
                                      Optimus.matrix_row_index,
                                      Optimus.matrix_col_index,
                                      Optimus.cell_calls,
                                      Optimus.loom_output_file,
  ])

  # Collect all of the pipeline metrics into a single Array
  Array[String] pipeline_metrics = [  Optimus.cell_metrics,
                                      Optimus.gene_metrics
  ]

  # Copy results of pipeline to test results bucket
  call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
      vault_token_path          = vault_token_path,
      google_account_vault_path = google_account_vault_path,
      destination_cloud_path    = results_path
  }

  # If updating truth then copy pipeline results to truth bucket
  if (update_truth){
    call Copy.CopyFilesFromCloudToCloud as CopyToTruth {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
      vault_token_path          = vault_token_path,
      google_account_vault_path = google_account_vault_path,
      destination_cloud_path    = truth_path
    }
  }

  # If not updating truth then gather the inputs and call verification wdl
  if (!update_truth) {
    call Utilities.GetValidationInputs as GetLoomInputs {
      input:
        input_file   = Optimus.loom_output_file,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetBamInputs {
      input:
        input_file   = Optimus.bam,
        results_path = results_path,
        truth_path   = truth_path,
    }

    call Utilities.GetValidationInputs as GetGeneMetrics {
      input:
        input_file   = Optimus.gene_metrics,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetCellMetrics {
      input:
        input_file   = Optimus.cell_metrics,
        results_path = results_path,
        truth_path   = truth_path
    }

    call VerifyOptimus.VerifyOptimus as Verify {
      input:
        test_loom          = GetLoomInputs.results_file,
        test_bam           = GetBamInputs.results_file,
        test_gene_metrics  = GetGeneMetrics.results_file,
        test_cell_metrics  = GetCellMetrics.results_file,
        truth_loom         = GetLoomInputs.truth_file,
        truth_bam          = GetBamInputs.truth_file,
        truth_gene_metrics = GetGeneMetrics.truth_file,
        truth_cell_metrics = GetCellMetrics.truth_file,
        done               = CopyToTestResults.done
    }

  }

  output {

  }

}