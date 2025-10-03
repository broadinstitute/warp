version 1.0


import "../../pipelines/wdl/dna_seq/somatic/single_sample/wgs/gdc_genome/GDCWholeGenomeSomaticSingleSample.wdl" as GDCWholeGenomeSomaticSingleSample
import "../../verification/VerifyGDCSomaticSingleSample.wdl" as VerifyGDCSomaticSingleSample
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestGDCWholeGenomeSomaticSingleSample {

    input {
      File? input_cram
      File? input_bam
      File? cram_ref_fasta
      File? cram_ref_fasta_index
      File? output_map
      String? unmapped_bam_suffix
      String base_file_name
      File? ubam
      File contamination_vcf
      File contamination_vcf_index
      File dbsnp_vcf
      File dbsnp_vcf_index
      File ref_fasta
      File ref_fai
      File ref_dict
      File ref_amb
      File ref_ann
      File ref_bwt
      File ref_pac
      File ref_sa

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      String vault_token_path
      String google_account_vault_path
    }

    meta {
      allowNestedInputs: true
    }
  
    call GDCWholeGenomeSomaticSingleSample.GDCWholeGenomeSomaticSingleSample {
      input:
        input_cram = input_cram,
        input_bam = input_bam,
        cram_ref_fasta = cram_ref_fasta, 
        cram_ref_fasta_index = cram_ref_fasta_index,
        output_map = output_map,
        unmapped_bam_suffix = unmapped_bam_suffix,
        base_file_name = base_file_name,
        ubam = ubam,
        contamination_vcf = contamination_vcf,
        contamination_vcf_index = contamination_vcf_index,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa
    }

    Array[String] pipeline_outputs = flatten([
                          [ # File outputs
                          GDCWholeGenomeSomaticSingleSample.bam,
                          GDCWholeGenomeSomaticSingleSample.bai,
                          GDCWholeGenomeSomaticSingleSample.insert_size_histogram_pdf,
                          GDCWholeGenomeSomaticSingleSample.contamination
                          ],
                          # Array[File]? outputs
                          flatten(select_all([
                            GDCWholeGenomeSomaticSingleSample.validation_report,
                            GDCWholeGenomeSomaticSingleSample.unmapped_bams
                          ]))
    ])

    Array[String] pipeline_metrics = flatten([
                          [ # File outputs
                          GDCWholeGenomeSomaticSingleSample.md_metrics,
                          GDCWholeGenomeSomaticSingleSample.insert_size_metrics
                          ]
    ])

    

    

    # Copy results of pipeline to test results bucket
    call Copy.CopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
        vault_token_path          = vault_token_path,
        google_account_vault_path = google_account_vault_path,
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.CopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
          vault_token_path          = vault_token_path,
          google_account_vault_path = google_account_vault_path,
          destination_cloud_path    = truth_path
      }
    }

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetMetrics {
          input:
            input_files = pipeline_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetBam {
          input:
            input_file = GDCWholeGenomeSomaticSingleSample.bam,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetBai {
          input:
            input_file = GDCWholeGenomeSomaticSingleSample.bai,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyGDCSomaticSingleSample.VerifyGDCSomaticSingleSample as Verify {
        input:
          truth_metrics = GetMetrics.truth_files, 
          test_metrics = GetMetrics.results_files,
          truth_bam = GetBam.truth_file, 
          test_bam = GetBam.results_file,
          truth_bai = GetBai.truth_file, 
          test_bai = GetBai.results_file,
          done = CopyToTestResults.done
      }
    }
}