version 1.0

import "../../../tasks/wdl/slidetags_utils/spatial-count.wdl" as SpatialCount
import "../../../tasks/wdl/slidetags_utils/positioning.wdl" as Positioning
import "../optimus/Optimus.wdl" as optimus

workflow SlideTags {

    String pipeline_version = "1.0.4"

    input {

        # slide-tags inputs
        Array[String] spatial_fastq
        Array[String] pucks
        # Dropsift is off by default
        Boolean run_dropsift = false

        # Optimus inputs
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

        String docker = "us.gcr.io/broad-gotc-prod/slide-tags:1.2.0"
     }
    
    parameter_meta {
        spatial_fastq: "Array of paths to spatial fastq files"
        pucks: "Array of paths to puck files"
        docker: "Docker image to use"
    }
    
    # Call the optimus workflow
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
            gex_expected_cells = expected_cells
    } 
    
    call SpatialCount.count as spatial_count {
        input:
            fastq_paths = spatial_fastq,
            pucks = pucks,
            docker = docker,
            input_id = input_id
     }

    call Positioning.generate_positioning as positioning {
        input:
            rna_paths = [Optimus.h5ad_output_file, Optimus.library_metrics],
            sb_path = spatial_count.sb_counts,
            docker = docker,
            input_id = input_id,
            run_dropsift = run_dropsift
     }

    output {
        # Version of Optimus pipeline
        String optimus_pipeline_version_out = Optimus.pipeline_version_out
        File optimus_genomic_reference_version = Optimus.genomic_reference_version
   
        # Optimus Metrics outputs
        File optimus_cell_metrics = Optimus.cell_metrics
        File optimus_gene_metrics = Optimus.gene_metrics
        File? optimus_cell_calls = Optimus.cell_calls
   
        # Optimus Star outputs 
        File optimus_library_metrics = Optimus.library_metrics
        File optimus_bam = Optimus.bam
        File optimus_matrix = Optimus.matrix
        File optimus_matrix_row_index = Optimus.matrix_row_index
        File optimus_matrix_col_index = Optimus.matrix_col_index
        File? optimus_aligner_metrics = Optimus.aligner_metrics
        File? optimus_mtx_files = Optimus.mtx_files
        File? optimus_filtered_mtx_files = Optimus.filtered_mtx_files
        File? optimus_multimappers_EM_matrix = Optimus.multimappers_EM_matrix
        File? optimus_multimappers_Uniform_matrix = Optimus.multimappers_Uniform_matrix
        File? optimus_multimappers_Rescue_matrix = Optimus.multimappers_Rescue_matrix
        File? optimus_multimappers_PropUnique_matrix = Optimus.multimappers_PropUnique_matrix
    
        # Optimus H5ad
        File optimus_h5ad_output_file = Optimus.h5ad_output_file
        
        # Optimus Cellbender outputs
        File? cb_cell_barcodes_csv = Optimus.cell_barcodes_csv
        File? cb_checkpoint_file = Optimus.checkpoint_file
        Array[File]? cb_h5_array = Optimus.h5_array
        Array[File]? cb_html_report_array = Optimus.html_report_array
        File? cb_log = Optimus.log
        Array[File]? cb_metrics_csv_array = Optimus.metrics_csv_array
        String? cb_output_directory = Optimus.output_directory
        File? cb_summary_pdf = Optimus.summary_pdf

        # Spatial/Positioning outputs
        File spatial_output_h5 = spatial_count.sb_counts
        File spatial_output_log = spatial_count.spatial_log
        File positioning_seurat_qs = positioning.seurat_qs
        File positioning_coords_csv = positioning.coords_csv
        File positioning_coords2_csv = positioning.coords2_csv
        File positioning_summary_pdf = positioning.summary_pdf
        File positioning_intermediates = positioning.intermediates_file
        File positioning_positioning_log = positioning.positioning_log
     }
}

