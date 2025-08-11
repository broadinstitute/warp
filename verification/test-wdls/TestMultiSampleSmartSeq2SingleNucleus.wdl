version 1.0


import "../../pipelines/wdl/smartseq2_single_nucleus_multisample/MultiSampleSmartSeq2SingleNucleus.wdl" as MultiSampleSmartSeq2SingleNucleus
import "../../verification/VerifyMultiSampleSmartSeq2SingleNucleus.wdl" as VerifyMultiSampleSmartSeq2SingleNucleus
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestMultiSampleSmartSeq2SingleNucleus {

    input {
      File genome_ref_fasta
      File tar_star_reference
      File annotations_gtf
      File adapter_list
      Array[String] input_ids
      Array[String]? input_names
      Array[String] fastq1_input_files
      Array[String] fastq2_input_files
      String batch_id
      String? batch_name
      Array[String]? project_id
      Array[String]? project_name
      Array[String]? library
      Array[String]? species
      Array[String]? organ
      String? input_name_metadata_field
      String? input_id_metadata_field

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth

      String cloud_provider
    }

    meta {
      allowNestedInputs: true
    }
  
    call MultiSampleSmartSeq2SingleNucleus.MultiSampleSmartSeq2SingleNucleus {
      input:
        genome_ref_fasta = genome_ref_fasta,
        tar_star_reference = tar_star_reference,
        annotations_gtf = annotations_gtf,
        adapter_list = adapter_list,
        input_ids = input_ids,
        input_names = input_names,
        fastq1_input_files = fastq1_input_files,
        fastq2_input_files = fastq2_input_files,
        batch_id = batch_id,
        batch_name = batch_name,
        project_id = project_id,
        project_name = project_name,
        library = library,
        species = species,
        organ = organ,
        input_name_metadata_field = input_name_metadata_field,
        input_id_metadata_field = input_id_metadata_field,
        cloud_provider = cloud_provider
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    MultiSampleSmartSeq2SingleNucleus.h5ad_output,
                                    ],
                                    # Array[File] outputs
                                    MultiSampleSmartSeq2SingleNucleus.bam_files,
                                    MultiSampleSmartSeq2SingleNucleus.exon_intron_count_files,
                                    
    ])

    

    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs]),
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs]),
          destination_cloud_path    = truth_path
      }
    }

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetBams {
          input:
            input_files = MultiSampleSmartSeq2SingleNucleus.bam_files,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetH5ad {
          input:
            input_file = MultiSampleSmartSeq2SingleNucleus.h5ad_output,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyMultiSampleSmartSeq2SingleNucleus.VerifyMultiSampleSmartSeq2SingleNucleus as Verify {
        input:
          truth_bams = GetBams.truth_files, 
          test_bams = GetBams.results_files,
          truth_h5ad = GetH5ad.truth_file,
          test_h5ad = GetH5ad.results_file,
          done = CopyToTestResults.done
      }
    }
}