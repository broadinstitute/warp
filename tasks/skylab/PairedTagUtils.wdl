version 1.0

task PairedTagDemultiplex {

  input {
    #runtime values
    #String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.6-1692962087"
    String docker = "us.gcr.io/broad-gotc-prod/warp-tools:2.0.0"
    # name of the sample
    String input_id
    # user provided id
    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field
    # gene annotation file in GTF format
    File annotation_file
    # the file "merged-cell-metrics.csv.gz" that contains the cellwise metrics
    File cell_metrics
    # the file "merged-gene-metrics.csv.gz" that contains the  genwise metrics
    File gene_metrics
    # file (.npz)  that contains the count matrix
    File sparse_count_matrix
    # file (.npy) that contains the array of cell barcodes
    File cell_id
    # file (.npy) that contains the array of gene names
    File gene_id
    # emptydrops output metadata
    File? empty_drops_result
    String counting_mode = "sc_rna"
    String add_emptydrops_data = "yes"


    String pipeline_version

    Int preemptible = 3
    Int disk = 200
    Int machine_mem_mb = 16000
    Int cpu = 4
  }

  meta {
    description: "This task will converts some of the outputs of Optimus pipeline into a h5ad file"
  }

  parameter_meta {
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -euo pipefail

    touch empty_drops_result.csv

    if [ "~{counting_mode}" == "sc_rna" ]; then
        python3 /warptools/scripts/create_h5ad_optimus.py \
          ~{if defined(empty_drops_result) then "--empty_drops_file  " + empty_drops_result  else "--empty_drops_file empty_drops_result.csv "  } \
          --add_emptydrops_data ~{add_emptydrops_data} \
          --annotation_file ~{annotation_file} \
          --cell_metrics ~{cell_metrics} \
          --gene_metrics ~{gene_metrics} \
          --cell_id ~{cell_id} \
          --gene_id  ~{gene_id} \
          --output_path_for_h5ad "~{input_id}" \
          --input_id ~{input_id} \
          ~{"--input_name " + input_name} \
          ~{"--input_id_metadata_field " + input_id_metadata_field} \
          ~{"--input_name_metadata_field " + input_name_metadata_field} \
          --count_matrix ~{sparse_count_matrix} \
          --expression_data_type "exonic" \
          --pipeline_version ~{pipeline_version}
    else
        python3 /warptools/scripts/create_snrna_optimus_full_h5ad.py \
          --annotation_file ~{annotation_file} \
          --cell_metrics ~{cell_metrics} \
          --gene_metrics ~{gene_metrics} \
          --cell_id ~{cell_id} \
          --gene_id  ~{gene_id} \
          --output_path_for_h5ad "~{input_id}" \
          --input_id ~{input_id} \
          ~{"--input_name " + input_name} \
          ~{"--input_id_metadata_field " + input_id_metadata_field} \
          ~{"--input_name_metadata_field " + input_name_metadata_field} \
          --count_matrix ~{sparse_count_matrix} \
          --expression_data_type "whole_transcript"\
          --pipeline_version ~{pipeline_version}
    fi
  >>>

  runtime {
    docker: docker
    cpu: cpu  # note that only 1 thread is supported by pseudobam
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    disk: disk + " GB" # TES
    preemptible: preemptible
  }

  output {
    File h5ad_output = "~{input_id}.h5ad"
  }
}
