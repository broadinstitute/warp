version 1.0

import "../../pipelines/skylab/atac/atac.wdl" as ATAC
import "../../verification/VerifyATAC.wdl" as VerifyATAC
import "../../tasks/broad/Utilities.wdl" as Utilities
import "../../tasks/broad/TerraCopyFilesFromCloudToCloud.wdl" as Copy

workflow TestATAC {

    input {
      # Fastq inputs
      Array[String] read1_fastq_gzipped
      Array[String] read2_fastq_gzipped
      Array[String] read3_fastq_gzipped

      # Output prefix/base name for all intermediate files and pipeline outputs
      String input_id
      String cloud_provider
      # Additional library aliquot ID
      String? atac_nhash_id

      #Expected cells from library preparation
      Int atac_expected_cells = 3000

      # Option for running files with preindex
      Boolean preindex = false
    
      # BWA ref
      File tar_bwa_reference
      # BWA machine type -- to select number of splits 
      Int num_threads_bwa = 128
      Int mem_size_bwa = 512
      String cpu_platform_bwa = "Intel Ice Lake"
      String vm_size

      # Text file containing chrom_sizes for genome build (i.e. hg38)
      File chrom_sizes
      #File for annotations for calculating ATAC TSSE
      File annotations_gtf
      # Whitelist
      File whitelist

      # TrimAdapters input
      String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
      String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"

      # These values will be determined and injected into the inputs by the scala test framework
      String truth_path
      String results_path
      Boolean update_truth
      Boolean run_cellbender = false

      # Option for running peak calling
      Boolean peak_calling = false

      # Optional aligned ATAC bam file
      File? aligned_ATAC_bam
    }

    meta {
      allowNestedInputs: true
    }
  
    call ATAC.ATAC {
      input:
        read1_fastq_gzipped = read1_fastq_gzipped,
        read2_fastq_gzipped = read2_fastq_gzipped,
        read3_fastq_gzipped = read3_fastq_gzipped,
        input_id = input_id,
        cloud_provider = cloud_provider,
        atac_nhash_id = atac_nhash_id,
        atac_expected_cells = atac_expected_cells,
        preindex = preindex,
        tar_bwa_reference = tar_bwa_reference,
        num_threads_bwa = num_threads_bwa,
        mem_size_bwa = mem_size_bwa,
        cpu_platform_bwa = cpu_platform_bwa,
        vm_size = vm_size,
        chrom_sizes = chrom_sizes,
        annotations_gtf = annotations_gtf,
        whitelist = whitelist,
        adapter_seq_read1 = adapter_seq_read1,
        adapter_seq_read3 = adapter_seq_read3,
        peak_calling = peak_calling,
        aligned_ATAC_bam = aligned_ATAC_bam
    }


    # Collect all of the pipeline outputs into single Array[String]
    Array[String] pipeline_outputs = flatten([
                                    [ # atac file outputs
                                    ATAC.fragment_file,
                                    ATAC.bam_aligned_output,
                                    ATAC.snap_metrics,
                                    ATAC.library_metrics_file
                                    ],                                  
    ])

    
    # Collect all of the pipeline metrics into single Array[String]
    Array[String] pipeline_metrics = flatten([
                                    [ # File outputs
                                    ATAC.fragment_file,
                                    ATAC.bam_aligned_output,
                                    ATAC.snap_metrics,
                                    ATAC.library_metrics_file
                                    ],
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
        call Utilities.GetValidationInputs as GetAtacBam {
          input:
            input_file = ATAC.bam_aligned_output,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetFragmentFile {
          input:
            input_file = ATAC.fragment_file,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetSnapMetrics {
          input:
            input_file = ATAC.snap_metrics,
            results_path = results_path,
            truth_path = truth_path
        }
        call Utilities.GetValidationInputs as GetAtacLibraryMetrics {
            input:
            input_file = ATAC.library_metrics_file,
            results_path = results_path,
            truth_path = truth_path
        }
        call VerifyATAC.VerifyATAC as Verify {
          input:
            truth_atac_bam = GetAtacBam.truth_file,
            test_atac_bam = GetAtacBam.results_file,
            truth_fragment_file = GetFragmentFile.truth_file,
            test_fragment_file = GetFragmentFile.results_file,
            truth_atac_h5ad = GetSnapMetrics.truth_file,
            test_atac_h5ad = GetSnapMetrics.results_file,
            truth_atac_library_metrics = GetAtacLibraryMetrics.truth_file,
            test_atac_library_metrics = GetAtacLibraryMetrics.results_file,
            done = CopyToTestResults.done
        }
    }
}
