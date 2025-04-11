version 1.0

import "scripts/spatial-count.wdl" as SpatialCount
import "scripts/positioning.wdl" as Positioning
import "../../../pipelines/skylab/optimus/Optimus.wdl" as optimus

workflow SlideTags {

    String pipeline_version = "1.0.0"

    input {
        # slide-tags inputs
        String id
        Array[String] fastq_paths
        Array[String] pucks
        Array[String] rna_paths
        String sb_path

        # Optimus Inputs
        String cloud_provider = "gcp"
        String input_id
        Int expected_cells = 3000 ## copied from Multiome ?
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
        String star_strand_mode = "Reverse"
        Boolean count_exons = false
        File gex_whitelist
        String? soloMultiMappers
        String? gex_nhash_id

        String docker = "us.gcr.io/broad-gotc-prod/slide-tags:1.1.0"
     }
    
    parameter_meta {
        fastq_paths: "Array of paths to spatial fastq files"
        pucks: "Array of paths to puck files"
        docker: "Docker image to use"
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
            gex_expected_cells = expected_cells
    } 
    

    call SpatialCount.count as spatial_count {
        input:
            fastq_paths = fastq_paths,
            pucks = pucks,
            docker = docker,
            input_id = input_id
     }

    call Positioning.generate_positioning as positioning {
        input:
            rna_paths = rna_paths,
            sb_path = spatial_count.sb_counts,
            docker = docker,
            input_id = input_id
     }

    output {
        # Optimus outputs
        File genomic_reference_version_gex = Optimus.genomic_reference_version
        File bam_gex = Optimus.bam
        File matrix_gex = Optimus.matrix
        File matrix_row_index_gex = Optimus.matrix_row_index
        File matrix_col_index_gex = Optimus.matrix_col_index
        File cell_metrics_gex = Optimus.cell_metrics
        File gene_metrics_gex = Optimus.gene_metrics
        File? cell_calls_gex = Optimus.cell_calls
        File h5ad_output_file_gex = JoinBarcodes.gex_h5ad_file
        File? multimappers_EM_matrix = Optimus.multimappers_EM_matrix
        File? multimappers_Uniform_matrix = Optimus.multimappers_Uniform_matrix
        File? multimappers_Rescue_matrix = Optimus.multimappers_Rescue_matrix
        File? multimappers_PropUnique_matrix = Optimus.multimappers_PropUnique_matrix
        File? gex_aligner_metrics = Optimus.aligner_metrics
        File? library_metrics = Optimus.library_metrics
        File? mtx_files = Optimus.mtx_files

         # Cellbender outputs
        File? cell_barcodes_csv = Optimus.cell_barcodes_csv
        File? checkpoint_file = Optimus.checkpoint_file
        Array[File]? h5_array = Optimus.h5_array
        Array[File]? html_report_array = Optimus.html_report_array
        File? log = Optimus.log
        Array[File]? metrics_csv_array = Optimus.metrics_csv_array
        String? output_directory = Optimus.output_directory
        File? summary_pdf = Optimus.summary_pdf

        # Spatial/Positioning outputs
        File spatial_output_h5 = spatial_count.sb_counts
        File spatial_output_log = spatial_count.spatial_log
        File positioning_seurat_qs = positioning.seurat_qs
        File positioning_coords_csv = positioning.coords_csv
        File positioning_summary_pdf = positioning.summary_pdf
        File positioning_output_file = positioning.output_file
        File positioning_positioning_log = positioning.positioning_log

     }
    
}

