version 1.0

import "../../../pipelines/skylab/multiome/atac.wdl" as atac
import "../../../pipelines/skylab/optimus/Optimus.wdl" as optimus
import "../../../tasks/skylab/H5adUtils.wdl" as H5adUtils
import "https://raw.githubusercontent.com/broadinstitute/CellBender/v0.3.1/wdl/cellbender_remove_background.wdl" as CellBender

workflow Multiome {
    String pipeline_version = "2.2.2"

    input {
        String input_id

        # Optimus Inputs
        String counting_mode = "sn_rna"
        Array[File] gex_r1_fastq
        Array[File] gex_r2_fastq
        Array[File]? gex_i1_fastq        
        File tar_star_reference
        File annotations_gtf
        File ref_genome_fasta
        File? mt_genes
        Int tenx_chemistry_version = 3
        Int emptydrops_lower = 100
        Boolean force_no_check = false
        Boolean ignore_r1_read_length = false
        String star_strand_mode = "Forward"
        Boolean count_exons = false
        File gex_whitelist = "gs://gcp-public-data--broad-references/RNA/resources/arc-v1/737K-arc-v1_gex.txt"

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

        # CellBender
        Boolean run_cellbender = false

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
            ref_genome_fasta = ref_genome_fasta,
            mt_genes = mt_genes,
            tenx_chemistry_version = tenx_chemistry_version,
            whitelist = gex_whitelist,
            emptydrops_lower = emptydrops_lower,
            force_no_check = force_no_check,
            ignore_r1_read_length = ignore_r1_read_length,
            star_strand_mode = star_strand_mode,
            count_exons = count_exons,
    }

    # Call the ATAC workflow
    call atac.ATAC as Atac {
        input:
            read1_fastq_gzipped = atac_r1_fastq,
            read2_fastq_gzipped = atac_r2_fastq,
            read3_fastq_gzipped = atac_r3_fastq,
            input_id = input_id + "_atac",
            tar_bwa_reference = tar_bwa_reference,
            annotations_gtf = annotations_gtf,
            chrom_sizes = chrom_sizes,
            whitelist = atac_whitelist,
            adapter_seq_read1 = adapter_seq_read1,
            adapter_seq_read3 = adapter_seq_read3
    }
    call H5adUtils.JoinMultiomeBarcodes as JoinBarcodes {
        input:
            atac_h5ad = Atac.snap_metrics,
            gex_h5ad = Optimus.h5ad_output_file,
            gex_whitelist = gex_whitelist,
            atac_whitelist = atac_whitelist
    }

    # Call CellBender
    if (run_cellbender) {
        call CellBender.run_cellbender_remove_background_gpu as CellBender {
            input:
                sample_name = sample_name,
                input_file_unfiltered = input_file_unfiltered,
                barcodes_file = barcodes_file,
                genes_file = genes_file,
                checkpoint_file = checkpoint_file,
                truth_file = truth_file,
                output_bucket_base_directory = output_bucket_base_directory,
                docker_image = "us.gcr.io/broad-dsde-methods/cellbender:0.3.0",
                dev_git_hash__ =  dev_git_hash__,
                expected_cells = expected_cells,
                total_droplets_included = total_droplets_included,
                force_cell_umi_prior = force_cell_umi_prior,
                force_empty_umi_prior = force_empty_umi_prior,
                model = model,
                low_count_threshold = low_count_threshold,
                fpr = fpr,
                epochs = epochs,
                z_dim = z_dim,
                z_layers = z_layers,
                empty_drop_training_fraction = empty_drop_training_fraction,
                learning_rate = learning_rate,
                exclude_feature_types = exclude_feature_types,
                ignore_features = ignore_features,
                projected_ambient_count_threshold = projected_ambient_count_threshold,
                checkpoint_mins = checkpoint_mins,
                final_elbo_fail_fraction = final_elbo_fail_fraction,
                epoch_elbo_fail_fraction = epoch_elbo_fail_fraction,
                num_training_tries = num_training_tries,
                learning_rate_retry_mult = learning_rate_retry_mult,
                posterior_batch_size = posterior_batch_size,
                estimator_multiple_cpu = estimator_multiple_cpu,
                constant_learning_rate = constant_learning_rate,
                debug = debug,
                hardware_zones = "us-east1-d us-east1-c us-central1-a us-central1-c us-west1-b",
                hardware_disk_size_GB = 50,
                hardware_boot_disk_size_GB = 20
        }
    }

    meta {
        allowNestedInputs: true
    }

    output {
        
        String multiome_pipeline_version_out = pipeline_version

        # atac outputs
        File bam_aligned_output_atac = Atac.bam_aligned_output
        File fragment_file_atac = Atac.fragment_file
        File snap_metrics_atac = JoinBarcodes.atac_h5ad_file

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
    }
}
