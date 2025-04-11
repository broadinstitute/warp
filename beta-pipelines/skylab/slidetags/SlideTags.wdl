version 1.0

import "scripts/spatial-count.wdl" as SpatialCount
import "scripts/positioning.wdl" as Positioning
import "../../../pipelines/skylab/optimus/Optimus.wdl" as optimus
import "../../../tasks/skylab/TarFiles.wdl" as TarFiles

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
    
    call TarFiles.tar_files as tar_files {
        input:
            task_outputs = [Optimus.genomic_reference_version, Optimus.bam, Optimus.matrix, Optimus.matrix_row_index, Optimus.matrix_col_index, Optimus.cell_metrics, Optimus.gene_metrics]
            optional_outputs = [Optimus.cell_calls, Optimus.multimappers_EM_matrix, Optimus.multimappers_Uniform_matrix, Optimus.multimappers_Rescue_matrix, Optimus.multimappers_PropUnique_matrix, Optimus.aligner_metrics, Optimus.library_metrics, Optimus.mtx_files, Optimus.cell_barcodes_csv, Optimus.checkpoint_file, Optimus.h5_array, Optimus.html_report_array, Optimus.log, Optimus.metrics_csv_array, Optimus.output_directory, Optimus.summary_pdf]
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
        File Optimus_output = tar_files.tarred_output
      
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

