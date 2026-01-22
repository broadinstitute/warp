version 1.0

import "../../../tasks/wdl/FastqProcessing.wdl" as FastqProcessing
import "../../../tasks/wdl/StarAlign.wdl" as StarAlign
import "../../../tasks/wdl/Metrics.wdl" as Metrics
import "../../../tasks/wdl/RunEmptyDrops.wdl" as RunEmptyDrops
import "../../../tasks/wdl/CheckInputs.wdl" as OptimusInputChecks
import "../../../tasks/wdl/H5adUtils.wdl" as H5adUtils
import "../../../tasks/wdl/Utilities.wdl" as utils
import "https://raw.githubusercontent.com/aawdeh/CellBender/aa-cbwithoutcuda/wdl/cellbender_remove_background_azure.wdl" as CellBender_no_cuda
import "https://raw.githubusercontent.com/broadinstitute/CellBender/v0.3.0/wdl/cellbender_remove_background.wdl" as CellBender

workflow Optimus {
  meta {
    description: "The optimus 3' pipeline processes 10x genomics sequencing data based on the v2 chemistry. It corrects cell barcodes and UMIs, aligns reads, marks duplicates, and returns data as alignments in BAM format and as counts in sparse matrix exchange format."
    allowNestedInputs: true
  }

  input {
    String cloud_provider

    # Mode for counting either "sc_rna" or "sn_rna"
    String counting_mode = "sc_rna"

    # Sequencing data inputs
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    String input_id
    # String for additional library aliquot ID
    String? gex_nhash_id
    String output_bam_basename = input_id
    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field
    # organism reference parameters
    File tar_star_reference
    File annotations_gtf
    File? mt_genes
    String? soloMultiMappers = "Uniform"
    Int gex_expected_cells = 3000

    # CellBender
    Boolean run_cellbender = false
    Int cellbender_memory_GB = 32

    # Chemistry options include: 2 or 3
    Int tenx_chemistry_version
    # Whitelist is selected based on the input tenx_chemistry_version
    File whitelist = checkOptimusInput.whitelist_out

    # read_structure is based on v2 or v3 chemistry
    String read_struct = checkOptimusInput.read_struct_out

    # Emptydrops lower cutoff
    Int emptydrops_lower = 100

    # Set to true to override input checks and allow pipeline to proceed with invalid input
    Boolean force_no_check = false
    
    # Check that tenx_chemistry_version matches the length of the read 1 fastq;
    # Set to true if you expect that r1_read_length does not match length of UMIs/barcodes for 10x chemistry v2 (26 bp) or v3 (28 bp).
    Boolean ignore_r1_read_length = false

    # Set to Forward, Reverse, or Unstranded to account for stranded library preparations (per STARsolo documentation)
    String star_strand_mode = "Forward"
    
    # Set to true to count reads aligned to exonic regions in sn_rna mode
    Boolean count_exons = false

    # Set starsolo disk size as adjustable parameter
    # Int disk_starsolo

    # this pipeline does not set any preemptible varibles and only relies on the task-level preemptible settings
    # you could override the tasklevel preemptible settings by passing it as one of the workflows inputs
    # for example: `"Optimus.StarAlign.preemptible": 3` will let the StarAlign task, which by default disables the
    # usage of preemptible machines, attempt to request for preemptible instance up to 3 times. 

  }

  # Version of this pipeline
  String pipeline_version = "8.0.4"

  # this is used to scatter matched [r1_fastq, r2_fastq, i1_fastq] arrays
  Array[Int] indices = range(length(r1_fastq))

  # 10x parameters
  File gcp_whitelist_v2 = "gs://gcp-public-data--broad-references/RNA/resources/737k-august-2016.txt"
  File gcp_whitelist_v3 = "gs://gcp-public-data--broad-references/RNA/resources/3M-febrary-2018.txt"
  File azure_whitelist_v2 = "https://datasetpublicbroadref.blob.core.windows.net/dataset/RNA/resources/737k-august-2016.txt?sv=2020-04-08&si=prod&sr=c&sig=DQxmjB4D1lAfOW9AxIWbXwZx6ksbwjlNkixw597JnvQ%3D"
  File azure_whitelist_v3 = "https://datasetpublicbroadref.blob.core.windows.net/dataset/RNA/resources/3M-febrary-2018.txt?sv=2020-04-08&si=prod&sr=c&sig=DQxmjB4D1lAfOW9AxIWbXwZx6ksbwjlNkixw597JnvQ%3D"

  # Takes the first read1 FASTQ from the inputs to check for chemistry match
  File r1_single_fastq = r1_fastq[0]

  # docker images
  String picard_cloud_docker = "picard-cloud:2.26.10"
  String pytools_docker = "pytools:1.0.0-1661263730"
  String empty_drops_docker = "empty-drops:1.0.1-4.2"
  String star_docker = "star:1.0.1-2.7.11a-1692706072"
  String warp_tools_docker = "warp-tools:2.6.1"
  String star_merge_docker = "star-merge-npz:1.3.0"
  String samtools_star = "samtools-star:1.0.0-1.11-2.7.11a-1731516196"
  String samtools_star_python = "samtools-star-python:1.0.0"

  #TODO how do we handle these?
  String alpine_docker = "alpine-bash@sha256:965a718a07c700a5204c77e391961edee37477634ce2f9cf652a8e4c2db858ff"
  String gcp_alpine_docker_prefix = "bashell/"
  String acr_alpine_docker_prefix = "dsppipelinedev.azurecr.io/"
  String alpine_docker_prefix = if cloud_provider == "gcp" then gcp_alpine_docker_prefix else acr_alpine_docker_prefix

  String ubuntu_docker = "ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
  String gcp_ubuntu_docker_prefix = "gcr.io/gcp-runtimes/"
  String acr_ubuntu_docker_prefix = "dsppipelinedev.azurecr.io/"
  String ubuntu_docker_prefix = if cloud_provider == "gcp" then gcp_ubuntu_docker_prefix else acr_ubuntu_docker_prefix

  String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
  String acr_docker_prefix = "dsppipelinedev.azurecr.io/"

  # choose docker prefix based on cloud provider
  String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

  # make sure either gcp or azr is supplied as cloud_provider
  if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
    call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
      input:
        message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
    }
  }

  parameter_meta {
    r1_fastq: "forward read, contains cell barcodes and molecule barcodes"
    r2_fastq: "reverse read, contains cDNA fragment generated from captured mRNA"
    i1_fastq: "(optional) index read, for demultiplexing of multiple samples on one flow cell."
    input_id: "name of sample matching this file, inserted into read group header"
    input_id_metadata_field: "String that describes the metadata field containing the input_id"
    input_name: "User provided sample name or cell_names"
    input_name_metadata_field: "String that describes the metadata field containing the input_name"
    tar_star_reference: "star genome reference"
    annotations_gtf: "gtf containing annotations for gene tagging (must match star reference)"
    whitelist: "10x genomics cell barcode allowlist"
    tenx_chemistry_version: "10X Genomics v2 (10 bp UMI) or v3 chemistry (12bp UMI)"
    force_no_check: "Set to true to override input checks and allow pipeline to proceed with invalid input"
    star_strand_mode: "STAR mode for handling stranded reads. Options are 'Forward', 'Reverse, or 'Unstranded.' Default is Forward."
  }

  call OptimusInputChecks.checkOptimusInput {
    input:
      force_no_check = force_no_check,
      counting_mode = counting_mode,
      count_exons = count_exons,
      gcp_whitelist_v2 = gcp_whitelist_v2,
      gcp_whitelist_v3 = gcp_whitelist_v3,
      azure_whitelist_v2 = azure_whitelist_v2,
      azure_whitelist_v3 = azure_whitelist_v3,
      tenx_chemistry_version = tenx_chemistry_version,
      r1_fastq = r1_single_fastq,
      ignore_r1_read_length = ignore_r1_read_length,
      cloud_provider = cloud_provider,
      alpine_docker_path = alpine_docker_prefix + alpine_docker
  }

  call StarAlign.STARGenomeRefVersion as ReferenceCheck {
    input:
      tar_star_reference = tar_star_reference,
      ubuntu_docker_path = ubuntu_docker_prefix + ubuntu_docker
  }

  call StarAlign.STARsoloFastq as STARsoloFastq {
      input:
        r1_fastq = r1_fastq,
        r2_fastq = r2_fastq,
        star_strand_mode = star_strand_mode,
        white_list = whitelist,
        tar_star_reference = tar_star_reference,
        chemistry = tenx_chemistry_version,
        counting_mode = counting_mode,
        count_exons = count_exons,
        input_id = input_id,
        output_bam_basename = output_bam_basename,
        soloMultiMappers = soloMultiMappers,
        samtools_star_docker_path = docker_prefix + samtools_star_python
    }
  
  call Metrics.CalculateGeneMetrics as GeneMetrics {
    input:
      bam_input = STARsoloFastq.bam_output,
      mt_genes = mt_genes,
      original_gtf = annotations_gtf,
      input_id = input_id,
      warp_tools_docker_path = docker_prefix + warp_tools_docker
  }

  call Metrics.CalculateCellMetrics as CellMetrics {
    input:
      bam_input = STARsoloFastq.bam_output,
      mt_genes = mt_genes,
      original_gtf = annotations_gtf,
      input_id = input_id,
      warp_tools_docker_path = docker_prefix + warp_tools_docker
  }

  if (counting_mode == "sc_rna"){
    call RunEmptyDrops.RunEmptyDrops {
      input:
        sparse_count_matrix = STARsoloFastq.sparse_counts,
        row_index = STARsoloFastq.row_index,
        col_index = STARsoloFastq.col_index,
        emptydrops_lower = emptydrops_lower,
        empty_drops_docker_path = docker_prefix + empty_drops_docker
    }
  }

  if (!count_exons) {
    call H5adUtils.OptimusH5adGeneration{
      input:
        input_id = input_id,
        gex_nhash_id = gex_nhash_id,
        expected_cells = gex_expected_cells,
        input_name = input_name,
        input_id_metadata_field = input_id_metadata_field,
        input_name_metadata_field = input_name_metadata_field,
        annotation_file = annotations_gtf,
        library_metrics = STARsoloFastq.library_metrics,
        cellbarcodes = STARsoloFastq.outputbarcodes,
        cell_metrics = CellMetrics.cell_metrics,
        gene_metrics = GeneMetrics.gene_metrics,
        sparse_count_matrix = STARsoloFastq.sparse_counts,
        cell_id = STARsoloFastq.row_index,
        gene_id = STARsoloFastq.col_index,
        empty_drops_result = RunEmptyDrops.empty_drops_result,
        counting_mode = counting_mode,
        pipeline_version = "Optimus_v~{pipeline_version}",
        warp_tools_docker_path = docker_prefix + warp_tools_docker
    }
  }
  if (count_exons  && counting_mode=="sn_rna") {
    call H5adUtils.SingleNucleusOptimusH5adOutput as OptimusH5adGenerationWithExons{
      input:
        input_id = input_id,
        gex_nhash_id = gex_nhash_id,
        expected_cells = gex_expected_cells,
        input_name = input_name,
        counting_mode = counting_mode,
        input_id_metadata_field = input_id_metadata_field,
        input_name_metadata_field = input_name_metadata_field,
        annotation_file = annotations_gtf,
        library_metrics = STARsoloFastq.library_metrics,
        cellbarcodes = STARsoloFastq.outputbarcodes,
        cell_metrics = CellMetrics.cell_metrics,
        gene_metrics = GeneMetrics.gene_metrics,
        sparse_count_matrix = STARsoloFastq.sparse_counts,
        cell_id = STARsoloFastq.row_index,
        gene_id = STARsoloFastq.col_index,
        sparse_count_matrix_exon = STARsoloFastq.sparse_counts,
        cell_id_exon = STARsoloFastq.row_index,
        gene_id_exon = STARsoloFastq.col_index,
        pipeline_version = "Optimus_v~{pipeline_version}",
        warp_tools_docker_path = docker_prefix + warp_tools_docker
    }
  }

  # Call CellBender
  if (run_cellbender) {
    if (cloud_provider == "gcp") {
      call CellBender.run_cellbender_remove_background_gpu as CellBender {
        input:
          sample_name = input_id,
          input_file_unfiltered = final_h5ad_output,
          hardware_boot_disk_size_GB = 20,
          hardware_cpu_count = 4,
          hardware_disk_size_GB = 50,
          hardware_gpu_type = "nvidia-tesla-t4",
          hardware_memory_GB = cellbender_memory_GB,
          hardware_preemptible_tries = 2,
          hardware_zones = "us-central1-a us-central1-c",
          nvidia_driver_version = "470.82.01"
      }
    }
    if (cloud_provider == "azure") {
      call CellBender_no_cuda.run_cellbender_remove_background_gpu as CellBender_no_cuda {
        input:
          sample_name = input_id,
          input_file_unfiltered = final_h5ad_output,
          hardware_boot_disk_size_GB = 20,
          hardware_cpu_count = 4,
          hardware_disk_size_GB = 50,
          hardware_gpu_type = "nvidia-tesla-t4",
          hardware_memory_GB = cellbender_memory_GB,
          hardware_preemptible_tries = 2,
          hardware_zones = "us-central1-a us-central1-c",
          nvidia_driver_version = "470.82.01"
      }
    }
  }

  File final_h5ad_output = select_first([OptimusH5adGenerationWithExons.h5ad_output, OptimusH5adGeneration.h5ad_output])
  File final_library_metrics = select_first([OptimusH5adGenerationWithExons.library_metrics, OptimusH5adGeneration.library_metrics])


  output {
    # version of this pipeline
    String pipeline_version_out = pipeline_version
    File genomic_reference_version = ReferenceCheck.genomic_ref_version
   
    # Metrics outputs
    File cell_metrics = CellMetrics.cell_metrics
    File gene_metrics = GeneMetrics.gene_metrics
    File? cell_calls = RunEmptyDrops.empty_drops_result
   
    # Star outputs 
    File library_metrics = final_library_metrics
    File bam = STARsoloFastq.bam_output
    File matrix = STARsoloFastq.sparse_counts
    File matrix_row_index = STARsoloFastq.row_index
    File matrix_col_index = STARsoloFastq.col_index
    File? aligner_metrics = STARsoloFastq.cell_reads_out
    File? mtx_files = STARsoloFastq.mtx_files
    File? filtered_mtx_files = STARsoloFastq.filtered_mtx_files
    File? multimappers_EM_matrix = STARsoloFastq.multimappers_EM_matrix
    File? multimappers_Uniform_matrix = STARsoloFastq.multimappers_Uniform_matrix
    File? multimappers_Rescue_matrix = STARsoloFastq.multimappers_Rescue_matrix
    File? multimappers_PropUnique_matrix = STARsoloFastq.multimappers_PropUnique_matrix
    
    # h5ad
    File h5ad_output_file = final_h5ad_output

    # cellbender outputs
    File? cell_barcodes_csv = CellBender.cell_csv
    File? checkpoint_file = CellBender.ckpt_file
    Array[File]? h5_array = CellBender.h5_array
    Array[File]? html_report_array = CellBender.report_array
    File? log = CellBender.log
    Array[File]? metrics_csv_array = CellBender.metrics_array
    String? output_directory = CellBender.output_dir
    File? summary_pdf = CellBender.pdf
  }
}