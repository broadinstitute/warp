version 1.0

import "../../../pipelines/skylab/atac/atac.wdl" as atac
import "../../../pipelines/skylab/optimus/Optimus.wdl" as optimus
import "../../../tasks/skylab/H5adUtils.wdl" as H5adUtils
import "../../../tasks/skylab/PairedTagUtils.wdl" as Demultiplexing
import "../../../tasks/broad/Utilities.wdl" as utils

workflow PairedTag {

    String pipeline_version = "1.8.3"


    input {
        String input_id
        # Additional library aliquot id
        String? gex_nhash_id
        String? atac_nhash_id

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
        File gex_whitelist = if cloud_provider == "gcp" then "gs://gcp-public-data--broad-references/RNA/resources/arc-v1/737K-arc-v1_gex.txt" else "https://datasetpublicbroadref.blob.core.windows.net/dataset/RNA/resources/arc-v1/737K-arc-v1_gex.txt?sv=2020-04-08&si=prod&sr=c&sig=DQxmjB4D1lAfOW9AxIWbXwZx6ksbwjlNkixw597JnvQ%3D"

        String? soloMultiMappers = "Uniform"
        # ATAC inputs
        # Array of input fastq files
        Array[File] atac_r1_fastq
        Array[File] atac_r2_fastq
        Array[File] atac_r3_fastq

        String vm_size = "Standard_M128s"

        # BWA input
        File tar_bwa_reference
        File chrom_sizes
        # Trimadapters input
        String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
        String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
        # Whitelist
        File atac_whitelist = if cloud_provider == "gcp" then "gs://gcp-public-data--broad-references/RNA/resources/arc-v1/737K-arc-v1_atac.txt" else "https://datasetpublicbroadref.blob.core.windows.net/dataset/RNA/resources/arc-v1/737K-arc-v1_atac.txt?sv=2020-04-08&si=prod&sr=c&sig=DQxmjB4D1lAfOW9AxIWbXwZx6ksbwjlNkixw597JnvQ%3D"

        # PairedTag
        Boolean preindex

        # Expected to be either 'gcp' or 'azure'
        String cloud_provider
    }

    # All docker images that are needed for tasks in this workflow
    String upstools_docker = "upstools:2.0.0"
    String snapatac_docker = "snapatac2:1.0.4-2.3.1-1700590229"

    # Prefixes based on cloud env
    String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
    String acr_docker_prefix = "dsppipelinedev.azurecr.io/"

    # choose docker prefix based on cloud_provider input
    String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

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
            cloud_provider = cloud_provider,
            soloMultiMappers = soloMultiMappers,
            gex_nhash_id = gex_nhash_id
    }

    # Call the ATAC workflow
        # Call the ATAC workflow
    scatter (idx in range(length(atac_r1_fastq))) {
        call Demultiplexing.PairedTagDemultiplex as demultiplex {
            input:
              read1_fastq = atac_r1_fastq[idx],
              read3_fastq = atac_r3_fastq[idx],
              barcodes_fastq = atac_r2_fastq[idx],
              input_id = input_id + "_atac",
              whitelist = atac_whitelist,
              preindex = preindex,
              docker_path = docker_prefix + upstools_docker
        }
    }

    # Call the ATAC workflow
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
            preindex = preindex,
            cloud_provider = cloud_provider,
            vm_size = vm_size,
            atac_nhash_id = atac_nhash_id
    }

    if (preindex) {
        call Demultiplexing.ParseBarcodes as ParseBarcodes {
            input:
              atac_h5ad = Atac_preindex.snap_metrics,
              atac_fragment = Atac_preindex.fragment_file,
              docker_path = docker_prefix + snapatac_docker,
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
        File atac_library_final = Atac_preindex.library_metrics_file

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
