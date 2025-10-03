version 1.0

import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../pipelines/wdl/optimus/Optimus.wdl" as Optimus
import "../../verification/VerifyOptimus.wdl" as VerifyOptimus
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestOptimus {

  input {

    # Mode for counting either "sc_rna" or "sn_rna"
    String counting_mode = "sc_rna"

    # Sequencing data inputs
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    String input_id
    String gex_nhash_id
    String output_bam_basename = input_id
    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field
    # organism reference parameters
    File tar_star_reference
    File annotations_gtf
    File? mt_genes
    String? soloMultiMappers

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

    # Set to Forward by default to count reads in 10x stranded mode
    String star_strand_mode

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

    String cloud_provider

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
      tenx_chemistry_version     = tenx_chemistry_version,
      emptydrops_lower           = emptydrops_lower,
      force_no_check             = force_no_check,
      star_strand_mode           = star_strand_mode,
      count_exons                = count_exons,
      ignore_r1_read_length      = ignore_r1_read_length,
      soloMultiMappers           = soloMultiMappers,
      cloud_provider             = cloud_provider,
      gex_nhash_id               = gex_nhash_id
  }

# Collect all of the pipeline outputs into single Array[String]
Array[String] pipeline_outputs = flatten([
                              [ # File outputs
                              Optimus.h5ad_output_file,
                              Optimus.matrix_col_index,
                              Optimus.matrix_row_index,
                              Optimus.matrix,
                              Optimus.bam,
                              Optimus.genomic_reference_version,
                              ],
                              # File? outputs
                              select_all([Optimus.mtx_files]),
                              select_all([Optimus.cell_calls]),
                              ])


  # Collect all of the pipeline metrics into single Array[String]
  Array[String] pipeline_metrics = flatten([
                              [ # File outputs
                              Optimus.gene_metrics,
                              Optimus.cell_metrics,
                              ],
                              # File? outputs
                              select_all([Optimus.library_metrics]),
                              select_all([Optimus.aligner_metrics]),
                              ])

  # Copy results of pipeline to test results bucket
  call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
      destination_cloud_path    = results_path
  }

  # If updating truth then copy pipeline results to truth bucket
  if (update_truth){
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
    input:
      files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
      destination_cloud_path    = truth_path
    }
  }

  # If not updating truth then gather the inputs and call verification wdl
  if (!update_truth) {
    call Utilities.GetValidationInputs as GetH5adInputs {
      input:
        input_file   = Optimus.h5ad_output_file,
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

  if(defined(Optimus.library_metrics)){
    call Utilities.GetValidationInputs as GetLibraryMetrics {
      input:
        input_file = Optimus.library_metrics,
        results_path = results_path,
        truth_path = truth_path
    }
}

    call VerifyOptimus.VerifyOptimus as Verify {
      input:
        test_h5ad          = GetH5adInputs.results_file,
        test_bam           = GetBamInputs.results_file,
        test_gene_metrics  = GetGeneMetrics.results_file,
        test_cell_metrics  = GetCellMetrics.results_file,
        truth_h5ad         = GetH5adInputs.truth_file,
        truth_bam          = GetBamInputs.truth_file,
        truth_gene_metrics = GetGeneMetrics.truth_file,
        truth_cell_metrics = GetCellMetrics.truth_file,
        test_library_metrics =  select_first([GetLibraryMetrics.results_file, ""]),
        truth_library_metrics = select_first([GetLibraryMetrics.truth_file, ""]),
        done               = CopyToTestResults.done
    }
  }

  output {

  }

}
