version 1.0


import "../../pipelines/wdl/paired_tag/PairedTag.wdl" as PairedTag
import "../../verification/VerifyPairedTag.wdl" as VerifyPairedTag
import "../../tasks/wdl/Utilities.wdl" as Utilities
import "../../tasks/wdl/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestPairedTag {

    input {
      String input_id
      String gex_nhash_id
      String atac_nhash_id

      # Optimus Inputs
      String counting_mode = "sn_rna"
      Array[File] gex_r1_fastq
      Array[File] gex_r2_fastq
      Array[File]? gex_i1_fastq
      File tar_star_reference
      File annotations_gtf
      File? mt_genes
      Int tenx_chemistry_version = 3
      Int emptydrops_lower = 100
      Boolean force_no_check = false
      Boolean ignore_r1_read_length = false
      String star_strand_mode = "Forward"
      Boolean count_exons = false
      File gex_whitelist = "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_gex.txt"
      String? soloMultiMappers = "Uniform"

      # ATAC inputs
      # Array of input fastq files
      Array[File] atac_r1_fastq
      Array[File] atac_r2_fastq
      Array[File] atac_r3_fastq

      Boolean preindex
      # BWA input
      File tar_bwa_reference

      # CreateFragmentFile input
      File chrom_sizes
      # Trimadapters input
      String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
      String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
      # Whitelist
      File atac_whitelist = "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_atac.txt"
      # Optional aligned ATAC bam file
      File? aligned_ATAC_bam

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      Boolean run_cellbender = false
      String cloud_provider

    }

    meta {
      allowNestedInputs: true
    }
  
    call PairedTag.PairedTag {
      input:
        counting_mode = counting_mode,
        gex_r1_fastq = gex_r1_fastq,
        gex_r2_fastq = gex_r2_fastq,
        gex_i1_fastq = gex_i1_fastq,
        input_id = input_id,
        tar_star_reference = tar_star_reference,
        annotations_gtf = annotations_gtf,
        mt_genes = mt_genes,
        tenx_chemistry_version = tenx_chemistry_version,
        emptydrops_lower = emptydrops_lower,
        force_no_check = force_no_check,
        ignore_r1_read_length = ignore_r1_read_length,
        star_strand_mode = star_strand_mode,
        count_exons = count_exons,
        gex_whitelist = gex_whitelist,
        preindex = preindex,
        atac_r1_fastq = atac_r1_fastq,
        atac_r2_fastq = atac_r2_fastq,
        atac_r3_fastq = atac_r3_fastq,
        tar_bwa_reference = tar_bwa_reference,
        adapter_seq_read1 = adapter_seq_read1,
        adapter_seq_read3 = adapter_seq_read3,
        chrom_sizes = chrom_sizes,
        atac_whitelist = atac_whitelist,
        soloMultiMappers = soloMultiMappers,
        cloud_provider = cloud_provider,
        gex_nhash_id = gex_nhash_id,
        atac_nhash_id = atac_nhash_id,
        aligned_ATAC_bam = aligned_ATAC_bam
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # Optimus file outputs
                                    PairedTag.h5ad_output_file_gex,
                                    PairedTag.matrix_col_index_gex,
                                    PairedTag.matrix_row_index_gex,
                                    PairedTag.matrix_gex,
                                    PairedTag.bam_gex,
                                    PairedTag.genomic_reference_version_gex,
                                    # atac file outputs
                                    PairedTag.fragment_file_atac,
                                    PairedTag.bam_aligned_output_atac,
                                    PairedTag.snap_metrics_atac
                                    ],
                                    # File? outputs
                                    select_all([PairedTag.cell_calls_gex]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    PairedTag.gene_metrics_gex,
                                    PairedTag.cell_metrics_gex,
                                    PairedTag.atac_library_final
                                    ],
                                    select_all([PairedTag.library_metrics]),
    ])

    # Copy results of pipeline to test results bucket
    call Copy.TerraCopyFilesFromCloudToCloud as CopyToTestResults {
      input:
        files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
        destination_cloud_path    = results_path
    }
  
    # If updating truth then copy output to truth bucket
    if (update_truth){
      call Copy.TerraCopyFilesFromCloudToCloud as CopyToTruth {
        input: 
          files_to_copy             = flatten([pipeline_outputs, pipeline_metrics]),
          destination_cloud_path    = truth_path
      }
    }

    # This is achieved by passing each desired file/array[files] to GetValidationInputs
    if (!update_truth){
        call Utilities.GetValidationInputs as GetOptimusH5ad {
          input:
            input_file = PairedTag.h5ad_output_file_gex,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetOptimusBam {
          input:
            input_file = PairedTag.bam_gex,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGeneMetrics {
          input:
            input_file = PairedTag.gene_metrics_gex,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCellMetrics {
          input:
            input_file = PairedTag.cell_metrics_gex,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetAtacBam {
          input:
            input_file = PairedTag.bam_aligned_output_atac,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFragmentFile {
          input:
            input_file = PairedTag.fragment_file_atac,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetSnapMetrics {
          input:
            input_file = PairedTag.snap_metrics_atac,
            results_path = results_path,
            truth_path = truth_path
        }

        if(defined(PairedTag.library_metrics)){
            call Utilities.GetValidationInputs as GetLibraryMetrics {
                input:
                    input_file = PairedTag.library_metrics,
                    results_path = results_path,
                    truth_path = truth_path
            }
        }

        call Utilities.GetValidationInputs as GetAtacLibraryMetrics {
            input:
                input_file = PairedTag.atac_library_final,
                results_path = results_path,
                truth_path = truth_path
        }

      call VerifyPairedTag.VerifyPairedTag as Verify {
        input:
          truth_optimus_h5ad = GetOptimusH5ad.truth_file,
          test_optimus_h5ad = GetOptimusH5ad.results_file,
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
          truth_atac_h5ad = GetSnapMetrics.truth_file,
          test_atac_h5ad = GetSnapMetrics.results_file,
          test_library_metrics =  select_first([GetLibraryMetrics.results_file, ""]),
          truth_library_metrics = select_first([GetLibraryMetrics.truth_file, ""]),
          test_atac_library_metrics = GetAtacLibraryMetrics.results_file,
          truth_atac_library_metrics = GetAtacLibraryMetrics.truth_file,
          done = CopyToTestResults.done
      }
    }
}