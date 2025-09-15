version 1.0

import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../pipelines/skylab/slidetags/SlideTags.wdl" as SlideTags
import "../../verification/VerifySlideTags.wdl" as VerifySlideTags
import "../../tasks/broad/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestSlideTags {

  input {

 	  #String id
      Array[String] spatial_fastq
      Array[String] pucks

      # Optimus Inputs
      Array[File] gex_r1_fastq
      Array[File] gex_r2_fastq
      Array[File]? gex_i1_fastq        
      File tar_star_reference
      File annotations_gtf
      File gex_whitelist
      String cloud_provider = "gcp"
      String input_id
      Int expected_cells = 3000 
      String counting_mode = "sn_rna"
      Int tenx_chemistry_version = 3
      Int emptydrops_lower = 100
      Boolean force_no_check = false
      Boolean ignore_r1_read_length = false
      String star_strand_mode = "Reverse"
      Boolean count_exons = false
      String? soloMultiMappers
      String? gex_nhash_id
	    File? mt_genes
	    String docker = "us.gcr.io/broad-gotc-prod/slide-tags@sha256:3f6d5e8612fbe32ca8cab09f32e1c83dec200a8382ca7cc23b48bb401d0212f4"
	
	    String truth_path
	    String results_path
	    Boolean update_truth

  }

  meta {
    allowNestedInputs: true
  }

  call SlideTags.SlideTags {
    input:
     # id              		    = id,
      spatial_fastq           = spatial_fastq,
      pucks                   = pucks,
      gex_r1_fastq            = gex_r1_fastq,
      gex_r2_fastq           	= gex_r2_fastq,
      gex_i1_fastq        	  = gex_i1_fastq,
      tar_star_reference      = tar_star_reference,
      annotations_gtf    	    = annotations_gtf,
      gex_whitelist  		      = gex_whitelist,
      cloud_provider         	= cloud_provider,
      input_id            	  = input_id,
      expected_cells     	    = expected_cells,
      counting_mode           = counting_mode,
      tenx_chemistry_version  = tenx_chemistry_version,
      emptydrops_lower        = emptydrops_lower,
      force_no_check          = force_no_check,
      ignore_r1_read_length   = ignore_r1_read_length,
      star_strand_mode        = star_strand_mode,
      count_exons             = count_exons,
      soloMultiMappers        = soloMultiMappers,
	    gex_nhash_id		        = gex_nhash_id,
	    mt_genes                = mt_genes,
	    docker                  = docker
  }

# collect all of the pipeline outputs into single Array[String]
Array[String] pipeline_outputs = flatten([
                              [ # File outputs
                              SlideTags.optimus_h5ad_output_file,
                              SlideTags.optimus_matrix_col_index,
                              SlideTags.optimus_matrix_row_index,
                              SlideTags.optimus_matrix,
                              SlideTags.optimus_bam,
                              SlideTags.optimus_genomic_reference_version,
					          SlideTags.spatial_output_h5,
					          SlideTags.positioning_seurat_qs,
					          SlideTags.positioning_coords_csv,
					          SlideTags.positioning_coords2_csv,
                              SlideTags.positioning_intermediates,
                              ],
                              # File? outputs
                              select_all([SlideTags.optimus_mtx_files]),
                              select_all([SlideTags.optimus_cell_calls]),
                              ])


  # Collect all of the pipeline metrics into single Array[String]
  Array[String] pipeline_metrics = flatten([
                              [ # File outputs
                              SlideTags.optimus_gene_metrics,
                              SlideTags.optimus_cell_metrics,
                              ],
                              # File? outputs
                              select_all([SlideTags.optimus_library_metrics]),
                              select_all([SlideTags.optimus_aligner_metrics]),
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
        input_file   = SlideTags.optimus_h5ad_output_file,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetBamInputs {
      input:
        input_file   = SlideTags.optimus_bam,
        results_path = results_path,
        truth_path   = truth_path,
    }

    call Utilities.GetValidationInputs as GetGeneMetrics {
      input:
        input_file   = SlideTags.optimus_gene_metrics,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetCellMetrics {
      input:
        input_file   = SlideTags.optimus_cell_metrics,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetSpatialOutput {
      input:
        input_file   = SlideTags.spatial_output_h5,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetSeuratQs {
      input:
        input_file   = SlideTags.positioning_seurat_qs,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetCoordsCsv {
      input:
        input_file   = SlideTags.positioning_coords_csv,
        results_path = results_path,
        truth_path   = truth_path
    }

    call Utilities.GetValidationInputs as GetCoords2Csv {
      input:
        input_file   = SlideTags.positioning_coords2_csv,
        results_path = results_path,
        truth_path   = truth_path
    }
    call Utilities.GetValidationInputs as GetIntermediates {
      input:
        input_file   = SlideTags.positioning_intermediates,
        results_path = results_path,
        truth_path   = truth_path
    }

  if(defined(SlideTags.optimus_library_metrics)){
    call Utilities.GetValidationInputs as GetLibraryMetrics {
      input:
        input_file = SlideTags.optimus_library_metrics,
        results_path = results_path,
        truth_path = truth_path
    }
}

    call VerifySlideTags.VerifySlideTags as Verify {
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
        test_spatial_output_h5 = GetSpatialOutput.results_file,
        truth_spatial_output_h5 = GetSpatialOutput.truth_file,
        test_seurat_qs = GetSeuratQs.results_file,
        truth_seurat_qs = GetSeuratQs.truth_file,
        test_coords_csv = GetCoordsCsv.results_file,
        truth_coords_csv = GetCoordsCsv.truth_file,
        test_coords2_csv = GetCoords2Csv.results_file,
        truth_coords2_csv = GetCoords2Csv.truth_file,
        test_intermediates_file = GetIntermediates.results_file,
        truth_intermediates_file = GetIntermediates.truth_file,
        done               = CopyToTestResults.done
    }
  }

  output {

  }

}
