version 1.0

import "../../beta-pipelines/multiome/atac.wdl" as atac
import "../../pipelines/skylab/optimus/Optimus.wdl" as optimus

workflow Multiome {
  String pipeline_version = "1.0.0"

  input {
      # Optimus Inputs
      String counting_mode = "sc_rna"
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
      String use_strand_info = "false"
      Boolean count_exons = false
      File whitelist= "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1.txt.gz"
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
    whitelist = whitelist,
    emptydrops_lower = emptydrops_lower,
    force_no_check = force_no_check,
    ignore_r1_read_length = ignore_r1_read_length,
    use_strand_info = use_strand_info,
    count_exons = count_exons
  }

  meta {
    allowNestedInputs: true
  }

  output {
    File bam = Optimus.bam
  }


}