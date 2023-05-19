version 1.0

import "../../../pipelines/skylab/multiome/atac.wdl" as atac
import "../../../pipelines/skylab/optimus/Optimus.wdl" as optimus

workflow Multiome {
    String pipeline_version = "1.0.0"

    input {
        # Optimus Inputs
        String counting_mode = "sn_rna"
        Array[File] r1_fastq
        Array[File] r2_fastq
        Array[File]? i1_fastq
        String input_id
        String output_bam_basename = input_id
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
        File gex_whitelist = "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_gex.txt"

        # ATAC inputs
        # Array of input fastq files
        Array[File] read1_fastq_gzipped
        Array[File] read2_fastq_gzipped
        Array[File] read3_fastq_gzipped
        # Output name
        String output_base_name
        # BWA input
        File tar_bwa_reference
        # CreateFragmentFile input
        Boolean barcodes_in_read_name
        File atac_gtf
        File chrom_sizes
        # Trimadapters input
        String adapter_seq_read1 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
        String adapter_seq_read3 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
        # Whitelist
        File atac_whitelist = "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_atac.txt"

        # script for monitoring tasks
        File monitoring_script

    }

    # Call the Optimus workflow
    call optimus.Optimus as Optimus {
        input:
            counting_mode = counting_mode,
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            i1_fastq = i1_fastq,
            input_id = input_id,
            output_bam_basename = output_bam_basename,
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
            count_exons = count_exons
    }

    # Call the ATAC workflow
    call atac.ATAC as Atac {
        input:
            read1_fastq_gzipped = read1_fastq_gzipped,
            read2_fastq_gzipped = read2_fastq_gzipped,
            read3_fastq_gzipped = read3_fastq_gzipped,
            output_base_name = output_base_name,
            tar_bwa_reference = tar_bwa_reference,
            monitoring_script = monitoring_script,
            barcodes_in_read_name = barcodes_in_read_name,
            atac_gtf = atac_gtf,
            chrom_sizes = chrom_sizes,
            whitelist = atac_whitelist,
            adapter_seq_read1 = adapter_seq_read1,
            adapter_seq_read3 = adapter_seq_read3
    }
    meta {
        allowNestedInputs: true
    }

    output {
        # atac outputs
        File bam_aligned_output = Atac.bam_aligned_output
        File fragment_file = Atac.fragment_file
        File snap_metrics = Atac.snap_metrics

        # optimus outputs
        String pipeline_version_out = Optimus.pipeline_version_out
        File genomic_reference_version = Optimus.genomic_reference_version
        File bam = Optimus.bam
        File matrix = Optimus.matrix
        File matrix_row_index = Optimus.matrix_row_index
        File matrix_col_index = Optimus.matrix_col_index
        File cell_metrics = Optimus.cell_metrics
        File gene_metrics = Optimus.gene_metrics
        File? cell_calls = Optimus.cell_calls
        File loom_output_file = Optimus.loom_output_file
    }
}