version 1.0


import "../../beta-pipelines/multiome/multiome.wdl" as Multiome
import "../../verification/VerifyMultiome.wdl" as VerifyMultiome
import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../tasks/broad/CopyFilesFromCloudToCloud.wdl" as Copy

workflow TestMultiome {

    input {
      String counting_mode = "sn_rna"
      Array[File] r1_fastq
      Array[File] r2_fastq
      Array[File]? i1_fastq
      String input_id
      String output_bam_basename = input_id
      File tar_star_reference
      File annotations_gtf
      File ref_genome_fasta
      File? mt_genes
      Int tenx_chemistry_version = 3
      Int emptydrops_lower = 100
      Boolean force_no_check = false
      Boolean ignore_r1_read_length = false
      String use_strand_info = "false"
      Boolean count_exons = false
      File gex_whitelist = "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1.txt"
      Array[File] read1_fastq_gzipped
      Array[File] read2_fastq_gzipped
      Array[File] read3_fastq_gzipped
      String output_base_name
      File tar_bwa_reference
      #File monitoring_script
      Boolean barcodes_in_read_name
      String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
      String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"

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
  
    call Multiome.Multiome {
      input:
        counting_mode = counting_mode,
        r1_fastq = r1_fastq,
        r2_fastq = r2_fastq,
        i1_fastq = i1_fastq,
        input_id = input_id,
        output_bam_basename = output_bam_basename,
        tar_star_reference = tar_star_reference,
        annotations_gtf = annotations_gtf,
        ref_genome_fasta = ref_genome_fasta,
        mt_genes = mt_genes,
        tenx_chemistry_version = tenx_chemistry_version,
        emptydrops_lower = emptydrops_lower,
        force_no_check = force_no_check,
        ignore_r1_read_length = ignore_r1_read_length,
        use_strand_info = use_strand_info,
        count_exons = count_exons,
        gex_whitelist = gex_whitelist,
        read1_fastq_gzipped = read1_fastq_gzipped,
        read2_fastq_gzipped = read2_fastq_gzipped,
        read3_fastq_gzipped = read3_fastq_gzipped,
        output_base_name = output_base_name,
        tar_bwa_reference = tar_bwa_reference,
        #monitoring_script = monitoring_script,
        barcodes_in_read_name = barcodes_in_read_name,
        adapter_seq_read1 = adapter_seq_read1,
        adapter_seq_read3 = adapter_seq_read3
  
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # File outputs
                                    Multiome.loom_output_file,
                                    Multiome.matrix_col_index,
                                    Multiome.matrix_row_index,
                                    Multiome.matrix,
                                    Multiome.bam,
                                    Multiome.genomic_reference_version,
                                    Multiome.fragment_file,
                                    Multiome.bam_aligned_output,
                                    ],
                                    # File? outputs
                                    select_all([Multiome.cell_calls]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    Multiome.gene_metrics,
                                    Multiome.cell_metrics,
                                    ],
                                    
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
        call Utilities.GetValidationInputs as GetLoom {
          input:
            input_file = Multiome.loom_output_file,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetOptimusBam {
          input:
            input_file = Multiome.bam,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGeneMetrics {
          input:
            input_file = Multiome.gene_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCellMetrics {
          input:
            input_file = Multiome.cell_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetAtacBam {
          input:
            input_file = Multiome.bam_aligned_output,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFragmentFile {
          input:
            input_file = Multiome.fragment_file,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyMultiome.VerifyMultiome as Verify {
        input:
          truth_loom = GetLoom.truth_file, 
          test_loom = GetLoom.results_file,
          truth_optimus_bam = GetOptimusBam.truth_file, 
          test_optimus_bam = GetOptimusBam.results_file,
          truth_gene_metrics = GetGeneMetrics.truth_file, 
          test_gene_metrics = GetGeneMetrics.results_file,
          truth_cell_metrics = GetCellMetrics.truth_file, 
          test_cell_metrics = GetCellMetrics.results_file,
          truth_atac_bam = GetAtacBam.truth_file, 
          test_atac_bam = GetAtacBam.results_file,
          truth_fragment_file = GetFragmentFile.truth_file,
          test_fragment_file = GetFragmentFile.results_file,
          done = CopyToTestResults.done
      }
    }
}