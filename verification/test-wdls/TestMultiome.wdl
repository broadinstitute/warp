version 1.0


import "../../pipelines/skylab/multiome/Multiome.wdl" as Multiome
import "../../verification/VerifyMultiome.wdl" as VerifyMultiome
import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../tasks/broad/TerraCopyFilesFromCloudToCloud.wdl" as Copy


workflow TestMultiome {

    input {
      String input_id
      String cloud_provider
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
      String? soloMultiMappers

      # ATAC inputs
      # Array of input fastq files
      Array[File] atac_r1_fastq
      Array[File] atac_r2_fastq
      Array[File] atac_r3_fastq

      # BWA input
      File tar_bwa_reference

      # CreateFragmentFile input
      File chrom_sizes
      # Trimadapters input
      String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
      String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"

      # Optional aligned ATAC bam file
      File? aligned_ATAC_bam

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      Boolean run_cellbender = false
      Boolean run_peak_calling = false

    }

    meta {
      allowNestedInputs: true
    }
  
    call Multiome.Multiome {
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
        atac_r1_fastq = atac_r1_fastq,
        atac_r2_fastq = atac_r2_fastq,
        atac_r3_fastq = atac_r3_fastq,
        tar_bwa_reference = tar_bwa_reference,
        adapter_seq_read1 = adapter_seq_read1,
        adapter_seq_read3 = adapter_seq_read3,
        chrom_sizes = chrom_sizes,
        run_cellbender = run_cellbender,
        soloMultiMappers = soloMultiMappers,
        cloud_provider = cloud_provider,
        gex_nhash_id = gex_nhash_id,
        atac_nhash_id = atac_nhash_id,
        run_peak_calling = run_peak_calling,
        aligned_ATAC_bam = aligned_ATAC_bam
    }

    
    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # Optimus file outputs
                                    Multiome.h5ad_output_file_gex,
                                    Multiome.matrix_col_index_gex,
                                    Multiome.matrix_row_index_gex,
                                    Multiome.matrix_gex,
                                    Multiome.bam_gex,
                                    Multiome.genomic_reference_version_gex,
                                    # atac file outputs
                                    Multiome.fragment_file_atac,
                                    Multiome.bam_aligned_output_atac,
                                    Multiome.snap_metrics_atac
                                    ],
                                    # File? outputs
                                    select_all([Multiome.cell_calls_gex]),
                                    
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    Multiome.gene_metrics_gex,
                                    Multiome.cell_metrics_gex,
                                    Multiome.atac_library_metrics
                                    ],
                                    # File? outputs
                                    select_all([Multiome.library_metrics]),
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
            input_file = Multiome.h5ad_output_file_gex,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetOptimusBam {
          input:
            input_file = Multiome.bam_gex,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetGeneMetrics {
          input:
            input_file = Multiome.gene_metrics_gex,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetCellMetrics {
          input:
            input_file = Multiome.cell_metrics_gex,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetAtacBam {
          input:
            input_file = Multiome.bam_aligned_output_atac,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFragmentFile {
          input:
            input_file = Multiome.fragment_file_atac,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetSnapMetrics {
          input:
            input_file = Multiome.snap_metrics_atac,
            results_path = results_path,
            truth_path = truth_path
        }

        if(defined(Multiome.library_metrics)){
            call Utilities.GetValidationInputs as GetLibraryMetrics {
                input:
                    input_file = Multiome.library_metrics,
                    results_path = results_path,
                    truth_path = truth_path
            }
        }

        call Utilities.GetValidationInputs as GetAtacLibraryMetrics {
            input:
            input_file = Multiome.atac_library_metrics,
            results_path = results_path,
            truth_path = truth_path
        }

      call VerifyMultiome.VerifyMultiome as Verify {
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