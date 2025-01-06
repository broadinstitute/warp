version 1.0

import "../../../pipelines/skylab/atac/atac.wdl" as atac
import "../../../pipelines/skylab/optimus/Optimus.wdl" as optimus
import "../../../tasks/skylab/H5adUtils.wdl" as H5adUtils
import "../../../tasks/broad/Utilities.wdl" as utils
import "https://raw.githubusercontent.com/aawdeh/CellBender/aa-cbwithoutcuda/wdl/cellbender_remove_background_azure.wdl" as CellBender_no_cuda
import "https://raw.githubusercontent.com/broadinstitute/CellBender/v0.3.0/wdl/cellbender_remove_background.wdl" as CellBender

workflow Multiome {

    String pipeline_version = "5.9.4"

    input {
        String cloud_provider
        String input_id
        # Additional library aliquot ID
        String? gex_nhash_id
        String? atac_nhash_id
        Int expected_cells = 3000

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
        # VM size used for several ATAC tasks
        String vm_size = "Standard_M128s"
        # BWA tar reference
        File tar_bwa_reference
        # Chromosone sizes 
        File chrom_sizes
        # Trimadapters input
        String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
        String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"

        # CellBender
        Boolean run_cellbender = false
    }

    # Determine docker prefix based on cloud provider
    String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
    String acr_docker_prefix = "dsppipelinedev.azurecr.io/"
    String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

    # Define docker images
    String snap_atac_docker_image = "snapatac2:1.0.4-2.3.1-1700590229"

    # Define all whitelist files
    File gcp_gex_whitelist = "gs://gcp-public-data--broad-references/RNA/resources/arc-v1/737K-arc-v1_gex.txt"
    File gcp_atac_whitelist = "gs://gcp-public-data--broad-references/RNA/resources/arc-v1/737K-arc-v1_atac.txt"
    File azure_gex_whitelist = "https://datasetpublicbroadref.blob.core.windows.net/dataset/RNA/resources/arc-v1/737K-arc-v1_gex.txt?sv=2020-04-08&si=prod&sr=c&sig=DQxmjB4D1lAfOW9AxIWbXwZx6ksbwjlNkixw597JnvQ%3D"
    File azure_atac_whitelist = "https://datasetpublicbroadref.blob.core.windows.net/dataset/RNA/resources/arc-v1/737K-arc-v1_atac.txt?sv=2020-04-08&si=prod&sr=c&sig=DQxmjB4D1lAfOW9AxIWbXwZx6ksbwjlNkixw597JnvQ%3D"

    # Determine which whitelist files to use based on cloud provider
    File gex_whitelist = if cloud_provider == "gcp" then gcp_gex_whitelist else azure_gex_whitelist
    File atac_whitelist = if cloud_provider == "gcp" then gcp_atac_whitelist else azure_atac_whitelist

    # Make sure either 'gcp' or 'azure' is supplied as cloud_provider input. If not, raise an error
    if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
        call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
            input:
                message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
        }
    }

    # Call the Optimus workflow
    call optimus.Optimus as Optimus {
        input:
            cloud_provider = cloud_provider,
            counting_mode = counting_mode,
            r1_fastq = gex_r1_fastq,
            r2_fastq = gex_r2_fastq,
            i1_fastq = gex_i1_fastq,
            input_id = input_id + "_gex",
            output_bam_basename = input_id + "_gex",
            gex_nhash_id = gex_nhash_id,
            tar_star_reference = tar_star_reference,
            annotations_gtf = annotations_gtf,
            mt_genes = mt_genes,
            tenx_chemistry_version = tenx_chemistry_version,
            whitelist = gex_whitelist,
            emptydrops_lower = emptydrops_lower,
            force_no_check = force_no_check,
            ignore_r1_read_length = ignore_r1_read_length,
            star_strand_mode = star_strand_mode,
            count_exons = count_exons,
            soloMultiMappers = soloMultiMappers,
            cloud_provider = cloud_provider,
            gex_expected_cells = expected_cells,
            run_cellbender = run_cellbender
    }

    # Call the ATAC workflow
    call atac.ATAC as Atac {
        input:
            cloud_provider = cloud_provider,
            read1_fastq_gzipped = atac_r1_fastq,
            read2_fastq_gzipped = atac_r2_fastq,
            read3_fastq_gzipped = atac_r3_fastq,
            input_id = input_id + "_atac",
            tar_bwa_reference = tar_bwa_reference,
            chrom_sizes = chrom_sizes,
            whitelist = atac_whitelist,
            adapter_seq_read1 = adapter_seq_read1,
            vm_size = vm_size,
            annotations_gtf = annotations_gtf,
            atac_nhash_id = atac_nhash_id,
            adapter_seq_read3 = adapter_seq_read3,
            atac_expected_cells = expected_cells
    }
    call H5adUtils.JoinMultiomeBarcodes as JoinBarcodes {
        input:
            docker_path = docker_prefix + snap_atac_docker_image,
            atac_h5ad = Atac.snap_metrics,
            gex_h5ad = Optimus.h5ad_output_file,
            gex_whitelist = gex_whitelist,
            atac_whitelist = atac_whitelist,
            atac_fragment = Atac.fragment_file
    }


    meta {
        allowNestedInputs: true
    }

    
    output {
        
        String multiome_pipeline_version_out = pipeline_version

        # atac outputs
        File bam_aligned_output_atac = Atac.bam_aligned_output
        File fragment_file_atac = JoinBarcodes.atac_fragment_tsv
        File fragment_file_index = JoinBarcodes.atac_fragment_tsv_index
        File snap_metrics_atac = JoinBarcodes.atac_h5ad_file
        File atac_library_metrics = Atac.library_metrics_file

        # optimus outputs
        File genomic_reference_version_gex = Optimus.genomic_reference_version
        File bam_gex = Optimus.bam
        File matrix_gex = Optimus.matrix
        File matrix_row_index_gex = Optimus.matrix_row_index
        File matrix_col_index_gex = Optimus.matrix_col_index
        File cell_metrics_gex = Optimus.cell_metrics
        File gene_metrics_gex = Optimus.gene_metrics
        File? cell_calls_gex = Optimus.cell_calls
        File h5ad_output_file_gex = JoinBarcodes.gex_h5ad_file
        Array[File?] multimappers_EM_matrix = Optimus.multimappers_EM_matrix
        Array[File?] multimappers_Uniform_matrix = Optimus.multimappers_Uniform_matrix
        Array[File?] multimappers_Rescue_matrix = Optimus.multimappers_Rescue_matrix
        Array[File?] multimappers_PropUnique_matrix = Optimus.multimappers_PropUnique_matrix
        File? gex_aligner_metrics = Optimus.aligner_metrics
        File? library_metrics = Optimus.library_metrics
        File? mtx_files = Optimus.mtx_files
        File? cell_barcodes_csv = Optimus.cell_barcodes_csv
        File? checkpoint_file = Optimus.checkpoint_file
        Array[File]? h5_array = Optimus.h5_array
        Array[File]? html_report_array = Optimus.html_report_array
        File? log = Optimus.log
        Array[File]? metrics_csv_array = Optimus.metrics_csv_array
        String? output_directory = Optimus.output_directory
        File? summary_pdf = Optimus.summary_pdf

    }
}
