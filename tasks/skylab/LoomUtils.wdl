version 1.0

task SmartSeq2LoomOutput {
  input {
    #runtime values
    String docker = "us.gcr.io/broad-gotc-prod/pytools:1.0.0-1661263730"
    # the gene count file "<input_id>_rsem.genes.results" in the task results folder call-RSEMExpression
    File rsem_gene_results
    # file named "<input_id>_QCs.csv" in the folder  "call-GroupQCOutputs/glob-*" of the the SS2  output
    Array[File] smartseq_qc_files
    # name of the sample
    String input_id
    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field

    String pipeline_version
    Int preemptible = 3
    Int disk = 200
    Int machine_mem_mb = 16000
    Int cpu = 4
  }

  meta {
    description: "This  task will converts some of the outputs of Smart Seq 2 pipeline into a loom file"
  }

  parameter_meta {
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -euo pipefail

    python3 /usr/gitc/create_loom_ss2.py \
       --qc_files ~{sep=' ' smartseq_qc_files} \
       --rsem_genes_results  ~{rsem_gene_results} \
       --output_loom_path  "~{input_id}.loom" \
       --input_id ~{input_id} \
       ~{"--input_name '" + input_name + "'"} \
       ~{"--input_id_metadata_field " + input_id_metadata_field} \
       ~{"--input_name_metadata_field " + input_name_metadata_field} \
       --pipeline_version ~{pipeline_version}
  }

  runtime {
    docker: docker
    cpu: cpu  # note that only 1 thread is supported by pseudobam
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    disk: disk + " GB" # TES
    preemptible: preemptible
  }

  output {
    File loom_output = "~{input_id}.loom"
  }
}


task OptimusLoomGeneration {

  input {
    #runtime values
    String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1681406657"
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
    description: "This task will converts some of the outputs of Optimus pipeline into a loom file"
  }

  parameter_meta {
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -euo pipefail

    touch empty_drops_result.csv

    if [ "~{counting_mode}" == "sc_rna" ]; then
        python3 /warptools/scripts/create_loom_optimus.py \
          ~{if defined(empty_drops_result) then "--empty_drops_file  " + empty_drops_result  else "--empty_drops_file empty_drops_result.csv "  } \
          --add_emptydrops_data ~{add_emptydrops_data} \
          --annotation_file ~{annotation_file} \
          --cell_metrics ~{cell_metrics} \
          --gene_metrics ~{gene_metrics} \
          --cell_id ~{cell_id} \
          --gene_id  ~{gene_id} \
          --output_path_for_loom "~{input_id}.loom" \
          --input_id ~{input_id} \
          ~{"--input_name " + input_name} \
          ~{"--input_id_metadata_field " + input_id_metadata_field} \
          ~{"--input_name_metadata_field " + input_name_metadata_field} \
          --count_matrix ~{sparse_count_matrix} \
          --expression_data_type "exonic" \
          --pipeline_version ~{pipeline_version}
    else
        python3 /warptools/scripts/create_snrna_optimus.py \
          --annotation_file ~{annotation_file} \
          --cell_metrics ~{cell_metrics} \
          --gene_metrics ~{gene_metrics} \
          --cell_id ~{cell_id} \
          --gene_id  ~{gene_id} \
          --output_path_for_loom "~{input_id}.loom" \
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
    File loom_output = "~{input_id}.loom"
  }
}


task AggregateSmartSeq2Loom {
    input {
        Array[File] loom_input
        String batch_id
        String? batch_name
        String? project_id
        String? project_name
        String? library
        String? species
        String? organ
        String pipeline_version
        String docker = "us.gcr.io/broad-gotc-prod/pytools:1.0.0-1661263730"
        Int disk = 200
        Int machine_mem_mb = 4000
        Int cpu = 1
    }

    meta {
      description: "aggregate the loom output"
    }

    command {
      set -e
      
      # Merge the loom files
      python3 /usr/gitc/ss2_loom_merge.py \
      --input-loom-files ~{sep=' ' loom_input} \
      --output-loom-file "~{batch_id}.loom" \
      --batch_id ~{batch_id} \
      ~{"--batch_name " + batch_name} \
      ~{"--project_id " + project_id} \
      ~{"--project_name " + project_name} \
      ~{"--library " + library} \
      ~{"--species " + species} \
      ~{"--organ " + organ} \
      --pipeline_version ~{pipeline_version}


    }

    output {
        File loom_output_file = "~{batch_id}.loom"
    }

    runtime {
      docker: docker
      cpu: cpu
      memory: "~{machine_mem_mb} MiB"
      disks: "local-disk ~{disk} HDD"
      disk: disk + " GB" # TES
      preemptible: 3
      maxRetries: 1
    }
}



task SingleNucleusOptimusLoomOutput {

    input {
        #runtime values
        String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1681406657"
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
        # the file "merged-gene-metrics.csv.gz" that contains the  genwise metrics
        File sparse_count_matrix_exon
        # file (.npy) that contains the array of cell barcodes
        File cell_id_exon
        # file (.npy) that contains the array of gene names
        File gene_id_exon

        String pipeline_version

        Int preemptible = 3
        Int disk = 200
        Int machine_mem_mb = 16000
        Int cpu = 4
    }

    meta {
        description: "This task will converts some of the outputs of Optimus pipeline into a loom file"
    }

    parameter_meta {
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command {
        set -euo pipefail

        python3 /warptools/scripts/create_snrna_optimus_counts.py \
        --annotation_file ~{annotation_file} \
        --cell_metrics ~{cell_metrics} \
        --gene_metrics ~{gene_metrics} \
        --count_matrix_1 ~{sparse_count_matrix} \
        --cell_id_1 ~{cell_id} \
        --gene_id_1  ~{gene_id} \
        --count_matrix_2 ~{sparse_count_matrix_exon} \
        --cell_id_2 ~{cell_id_exon} \
        --gene_id_2  ~{gene_id_exon} \
        --output_path_for_loom "~{input_id}.loom" \
        --input_id ~{input_id} \
        ~{"--input_name " + input_name} \
        ~{"--input_id_metadata_field " + input_id_metadata_field} \
        ~{"--input_name_metadata_field " + input_name_metadata_field} \
        --expression_data_type "whole_transcript" \
        --pipeline_version ~{pipeline_version}
    }

    runtime {
        docker: docker
        cpu: cpu  # note that only 1 thread is supported by pseudobam
        memory: "~{machine_mem_mb} MiB"
        disks: "local-disk ~{disk} HDD"
        disk: disk + " GB" # TES
        preemptible: preemptible
    }

    output {
        File loom_output = "~{input_id}.loom"
    }
}



task SlideSeqLoomOutput {
  input {
    File bead_locations
    File cell_metrics
    File gene_metrics
    File cell_id
    File gene_id
    File sparse_count_matrix
    File annotation_file
    String input_id
    String pipeline_version

    String docker = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1681406657"
    Int disk_size_gb = 200
    Int memory_mb = 18000
    Int cpu = 4
  }

  command <<<
    python3 /warptools/scripts/create_loom_slide_seq.py \
       --bead_locations ~{bead_locations} \
       --annotation_file ~{annotation_file} \
       --cell_metrics ~{cell_metrics} \
       --gene_metrics ~{gene_metrics} \
       --cell_id ~{cell_id} \
       --gene_id  ~{gene_id} \
       --output_path_for_loom "~{input_id}.loom" \
       --input_id ~{input_id} \
       --count_matrix ~{sparse_count_matrix} \
       --pipeline_version ~{pipeline_version}
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File loom_output = "~{input_id}.loom"  }
}
