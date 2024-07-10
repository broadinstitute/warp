version 1.0

import "../../../pipelines/skylab/multiome/atac.wdl" as atac
import "../../../pipelines/skylab/optimus/Optimus.wdl" as optimus
import "../../../tasks/skylab/H5adUtils.wdl" as H5adUtils
import "../../../tasks/skylab/PairedTagUtils.wdl" as Demultiplexing
workflow PairedTag {
    String pipeline_version = "1.1.2"

    input {
        String input_id

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
        File gex_whitelist = "gs://gcp-public-data--broad-references/RNA/resources/arc-v1/737K-arc-v1_gex.txt"

        String? soloMultiMappers = "Uniform"
        # ATAC inputs
        # Array of input fastq files
        Array[File] atac_r1_fastq
        Array[File] atac_r2_fastq
        Array[File] atac_r3_fastq
        # BWA input
        File tar_bwa_reference
        File chrom_sizes
        # Trimadapters input
        String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
        String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
        # Whitelist
        File atac_whitelist = "gs://gcp-public-data--broad-references/RNA/resources/arc-v1/737K-arc-v1_atac.txt"

        # PairedTag
        Boolean preindex
    }
    # Call the Optimus workflow
    call optimus.Optimus as Optimus {
        input:
            counting_mode = counting_mode,
            r1_fastq = gex_r1_fastq,
            r2_fastq = gex_r2_fastq,
            i1_fastq = gex_i1_fastq,
            input_id = input_id + "_gex",
            output_bam_basename = input_id + "_gex",
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
            soloMultiMappers = soloMultiMappers
    }

    # Call the ATAC workflow
        # Call the ATAC workflow
    scatter (idx in range(length(atac_r1_fastq))) {
        call Demultiplexing.PairedTagDemultiplex as demultiplex {
            input:
              read1_fastq = atac_r1_fastq[idx],
              read3_fastq = atac_r3_fastq[idx],
              barcodes_fastq = atac_r2_fastq[idx],
              input_id = input_id,
              whitelist = atac_whitelist,
              preindex = preindex
        }
    }      
    call atac.ATAC as Atac_preindex {
        input:
            read1_fastq_gzipped = demultiplex.fastq1,
            read2_fastq_gzipped = demultiplex.barcodes,
            read3_fastq_gzipped = demultiplex.fastq3,
            input_id = input_id + "_atac",
            tar_bwa_reference = tar_bwa_reference,
            chrom_sizes = chrom_sizes,
            whitelist = atac_whitelist,
            adapter_seq_read1 = adapter_seq_read1,
            adapter_seq_read3 = adapter_seq_read3,
            annotations_gtf = annotations_gtf,
            preindex = preindex
    }

    if (preindex) {
        call Demultiplexing.ParseBarcodes as ParseBarcodes {
            input:
              atac_h5ad = Atac_preindex.snap_metrics,
              atac_fragment = Atac_preindex.fragment_file
        }
    }      

    meta {
        allowNestedInputs: true
    }
    
    File atac_fragment_out = select_first([ParseBarcodes.atac_fragment_tsv,Atac_preindex.fragment_file])
    File atac_h5ad_out = select_first([ParseBarcodes.atac_h5ad_file, Atac_preindex.snap_metrics])
    output {
        
        String pairedtag_pipeline_version_out = pipeline_version

        # atac outputs
        File bam_aligned_output_atac = Atac_preindex.bam_aligned_output
        File fragment_file_atac = atac_fragment_out
        File snap_metrics_atac = atac_h5ad_out

        # optimus outputs
        File genomic_reference_version_gex = Optimus.genomic_reference_version
        File bam_gex = Optimus.bam
        File matrix_gex = Optimus.matrix
        File matrix_row_index_gex = Optimus.matrix_row_index
        File matrix_col_index_gex = Optimus.matrix_col_index
        File cell_metrics_gex = Optimus.cell_metrics
        File gene_metrics_gex = Optimus.gene_metrics
        File? cell_calls_gex = Optimus.cell_calls
        File h5ad_output_file_gex = Optimus.h5ad_output_file
        File? library_metrics = Optimus.library_metrics
        Array[File?] multimappers_EM_matrix = Optimus.multimappers_EM_matrix
        Array[File?] multimappers_Uniform_matrix = Optimus.multimappers_Uniform_matrix
        Array[File?] multimappers_Rescue_matrix = Optimus.multimappers_Rescue_matrix
        Array[File?] multimappers_PropUnique_matrix = Optimus.multimappers_PropUnique_matrix
    }
}
