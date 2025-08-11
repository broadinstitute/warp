optimus_inputs = {
  "Optimus.annotations_gtf": "workspace.human_annotations_gtf",
  "Optimus.chemistry": "\"tenX_v2\"",
  "Optimus.counting_mode": "\"sc_rna\"",
  "Optimus.emptydrops_lower": "1",
  "Optimus.i1_fastq": "this.participant_lanes.i1_fastq",
  "Optimus.input_id": "this.participant_lane_set_id",
  "Optimus.input_id_metadata_field": "\"sequencing_process.provenance.document_id\"",
  "Optimus.input_name": "this.input_name",
  "Optimus.input_name_metadata_field": "\"sequencing_input.biomaterial_core.biomaterial_id\"",
  "Optimus.r1_fastq": "this.participant_lanes.r1_fastq",
  "Optimus.r2_fastq": "this.participant_lanes.r2_fastq",
  "Optimus.ref_genome_fasta": "workspace.human_ref_genome_fasta",
  "Optimus.tar_star_reference": "workspace.human_tar_star_reference",
  "Optimus.whitelist": "workspace.whitelist_v2"
  }


optimus_outputs = {
  "Optimus.bam": "this.bam",
  "Optimus.cell_calls": "this.cell_calls",
  "Optimus.cell_metrics": "this.cell_metrics",
  "Optimus.gene_metrics": "this.gene_metrics",
  "Optimus.loom_output_file": "this.loom_output_file",
  "Optimus.matrix": "this.matrix",
  "Optimus.matrix_col_index": "this.matrix_col_index",
  "Optimus.matrix_row_index": "this.matrix_row_index",
  "Optimus.pipeline_version_out": "this.pipeline_version_out"
  }

