version 1.0



task OptimusH5adGeneration {

  input {
    #runtime values
    String warp_tools_docker_path
    # name of the sample
    String input_id
    String? gex_nhash_id
    # user provided id
    String counting_mode = "sc_rna"
    Int expected_cells = 3000
    String? input_name
    String? input_id_metadata_field
    String? input_name_metadata_field
    # gene annotation file in GTF format
    File annotation_file
    File? cellbarcodes
    File? library_metrics
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
    #String counting_mode = "sc_rna"
    String add_emptydrops_data = "yes"
    String gtf_path = annotation_file


    String pipeline_version

    Int preemptible = 3
    Int disk = 200
    Int machine_mem_mb = 32000
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

    touch empty_drops_result.csvs

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
          --pipeline_version ~{pipeline_version} \
          --gtf_path ~{gtf_path}
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
          --pipeline_version ~{pipeline_version} \
          --gtf_path ~{gtf_path}
    fi

    # modify h5ad to include doublets, NHASHID, and build library metrics
    python3 /warptools/scripts/add_library_tso_doublets.py \
     --gex_h5ad "~{input_id}.h5ad" \
     --cellbarcodes ~{cellbarcodes} \
     ~{"--gex_nhash_id " + gex_nhash_id} \
     --library_csv ~{library_metrics} \
     --input_id ~{input_id} \
     --counting_mode ~{counting_mode} \
     --expected_cells ~{expected_cells}

    mv library_metrics.csv ~{input_id}_~{gex_nhash_id}_library_metrics.csv

  >>>

  runtime {
    docker: warp_tools_docker_path
    cpu: cpu  # note that only 1 thread is supported by pseudobam
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    disk: disk + " GB" # TES
    preemptible: preemptible
  }

  output {
    File h5ad_output = "~{input_id}.h5ad"
    File library_metrics = "~{input_id}_~{gex_nhash_id}_library_metrics.csv"
  }
}


task SingleNucleusOptimusH5adOutput {

    input {
        #runtime values
        String warp_tools_docker_path
        # name of the sample
        String input_id
        # additional aliquot id
        String? gex_nhash_id
        # user provided id
        String? counting_mode
        Int expected_cells = 3000
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
        # library-level metrics
        File? library_metrics
        # Cell calls from starsolo in TSV format
        File? cellbarcodes
        String gtf_path = annotation_file

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

        python3 /warptools/scripts/create_snrna_optimus_exons_h5ad.py \
        --annotation_file ~{annotation_file} \
        --cell_metrics ~{cell_metrics} \
        --gene_metrics ~{gene_metrics} \
        --count_matrix_1 ~{sparse_count_matrix} \
        --cell_id_1 ~{cell_id} \
        --gene_id_1  ~{gene_id} \
        --count_matrix_2 ~{sparse_count_matrix_exon} \
        --cell_id_2 ~{cell_id_exon} \
        --gene_id_2  ~{gene_id_exon} \
        --output_path_for_h5ad "~{input_id}" \
        --input_id ~{input_id} \
        ~{"--input_name " + input_name} \
        ~{"--input_id_metadata_field " + input_id_metadata_field} \
        ~{"--input_name_metadata_field " + input_name_metadata_field} \
        --expression_data_type "whole_transcript" \
        --pipeline_version ~{pipeline_version} \
        --gtf_path ~{gtf_path}

        # modify h5ad to include doublets, NHASHID, and build library metrics
        python3 /warptools/scripts/add_library_tso_doublets.py \
        --gex_h5ad "~{input_id}.h5ad" \
        --cellbarcodes ~{cellbarcodes} \
        ~{"--gex_nhash_id " + gex_nhash_id} \
        --library_csv ~{library_metrics} \
        --input_id ~{input_id} \
        --counting_mode ~{counting_mode} \
        --expected_cells ~{expected_cells}


        mv library_metrics.csv ~{input_id}_~{gex_nhash_id}_library_metrics.csv

    >>>
    runtime {
        docker: warp_tools_docker_path
        cpu: cpu  # note that only 1 thread is supported by pseudobam
        memory: "~{machine_mem_mb} MiB"
        disks: "local-disk ~{disk} HDD"
        disk: disk + " GB" # TES
        preemptible: preemptible
    }

    output {
        File h5ad_output = "~{input_id}.h5ad"
        File library_metrics = "~{input_id}_~{gex_nhash_id}_library_metrics.csv"
    }
}

task JoinMultiomeBarcodes {
  input {
    File atac_h5ad
    File atac_fragment
    File gex_h5ad
    File gex_whitelist
    File atac_whitelist
    String input_gtf
    String input_bwa_reference

    Int nthreads = 1
    String cpuPlatform = "Intel Cascade Lake"
    Int machine_mem_mb = ceil((size(atac_h5ad, "MiB") + size(gex_h5ad, "MiB") + size(atac_fragment, "MiB")) * 8) + 10000
    Int disk =  ceil((size(atac_h5ad, "GiB") + size(gex_h5ad, "GiB") + size(atac_fragment, "GiB")) * 8) + 100
    String docker_path
  }
  String gex_base_name = basename(gex_h5ad, ".h5ad")
  String atac_base_name = basename(atac_h5ad, ".h5ad")
  String atac_fragment_base = basename(atac_fragment, ".sorted.tsv.gz")

  parameter_meta {
    atac_h5ad: "The resulting h5ad from the ATAC workflow."
    atac_fragment: "The resulting fragment TSV from the ATAC workflow."
    gex_h5ad: "The resulting h5ad from the Optimus workflow."
    gex_whitelist: "Whitelist used for gene expression barcodes."
    atac_whitelist: "Whitelist used for ATAC barcodes."
    input_gtf: "Reference GTF file used in the analysis."
    input_bwa_reference: "Reference genome used in the analysis."
  }

  command <<<
    set -e pipefail

    # decompress the bgzipped fragment file
    echo "Moving fragment file for bgzipping"
    mv ~{atac_fragment} ~{atac_fragment_base}.sorted.tsv.gz
    echo "Decompressing fragment file"
    bgzip -d "~{atac_fragment_base}.sorted.tsv.gz"
    echo "Done decompressing"


    python3 <<CODE

    # set parameters
    atac_h5ad = "~{atac_h5ad}"
    atac_fragment = "~{atac_fragment_base}.sorted.tsv"
    gex_h5ad = "~{gex_h5ad}"
    gex_whitelist = "~{gex_whitelist}"
    atac_whitelist = "~{atac_whitelist}"
    input_gtf = "~{input_gtf}"
    input_bwa_reference = "~{input_bwa_reference}"

    # import anndata to manipulate h5ad files
    import anndata as ad
    import pandas as pd
    import snapatac2 as snap
    print("Reading ATAC h5ad:")
    print("~{atac_h5ad}")
    print("Read ATAC fragment file:")
    print(atac_fragment)
    print("Reading Optimus h5ad:")
    print("~{gex_h5ad}")
    atac_data = ad.read_h5ad("~{atac_h5ad}")
    gex_data = ad.read_h5ad("~{gex_h5ad}")
    atac_tsv = pd.read_csv(atac_fragment, sep="\t", names=['chr','start', 'stop', 'barcode','n_reads'])
    print("Printing ATAC fragment tsv")
    print(atac_tsv)
    whitelist_gex = pd.read_csv("~{gex_whitelist}", header=None, names=["gex_barcodes"])
    whitelist_atac = pd.read_csv("~{atac_whitelist}", header=None, names=["atac_barcodes"])

    # get dataframes
    df_atac = atac_data.obs
    df_gex = gex_data.obs
    print(df_atac)
    print(df_gex)

    # Idenitfy the barcodes in the whitelist that match barcodes in datasets
    print("Printing whitelist_gex")
    print(whitelist_gex[1:10])

    df_all = pd.concat([whitelist_gex,whitelist_atac], axis=1)
    df_both_gex = df_all.copy()
    df_both_atac = df_all.copy()
    df_both_atac.set_index("atac_barcodes", inplace=True)
    df_both_gex.set_index("gex_barcodes", inplace=True)
    df_atac = atac_data.obs.join(df_both_atac)
    df_gex = gex_data.obs.join(df_both_gex)
    df_fragment = pd.merge(atac_tsv, df_both_atac, left_on='barcode', right_index=True, how='left')
    # set atac_data.obs to new dataframe
    print("Setting ATAC obs to new dataframe")
    atac_data.obs = df_atac
    #rename ATAC matrix 'index' to atac_barcodes
    atac_data.obs.index.name = 'atac_barcodes'
    # set gene_data.obs to new dataframe
    print("Setting Optimus obs to new dataframe")
    gex_data.obs = df_gex

    import os

    # Add whitelist provenance metadata
    gex_data.uns["whitelists"] = {
    "gex_whitelist_gs_path": gex_whitelist,
    "atac_whitelist_gs_path": atac_whitelist
    }

    atac_data.uns["whitelists"] = {
    "gex_whitelist_gs_path": gex_whitelist,
    "atac_whitelist_gs_path": atac_whitelist
    }

    # write out the names of the whitelists in separate text files for provenance tracking
    gex_whitelist_name = os.path.basename(gex_whitelist)
    atac_whitelist_name = os.path.basename(atac_whitelist)

    with open("gex_whitelist_used.txt", "w") as f:
    f.write(gex_whitelist_name)

    with open("atac_whitelist_used.txt", "w") as f:
    f.write(atac_whitelist_name)

    # write out the files
    gex_data.write("~{gex_base_name}.h5ad")
    atac_data.write_h5ad("~{atac_base_name}.h5ad")
    df_fragment.to_csv("~{atac_fragment_base}.tsv", sep='\t', index=False, header = False)
    CODE

    # Add reference information to fragment file header
    echo "Adding reference information to fragment file"
    echo "# Reference genome is ~{input_bwa_reference}" > "~{atac_fragment_base}.with_header.tsv"
    echo "# Reference GTF is ~{input_gtf}" >> "~{atac_fragment_base}.with_header.tsv"
    cat "~{atac_fragment_base}.tsv" >> "~{atac_fragment_base}.with_header.tsv"
    mv "~{atac_fragment_base}.with_header.tsv" "~{atac_fragment_base}.tsv"

    # sorting the file (skip header lines that start with #)
    echo "Sorting file"
    (head -n 2 "~{atac_fragment_base}.tsv"; tail -n +3 "~{atac_fragment_base}.tsv" | sort -k1,1V -k2,2n) > "~{atac_fragment_base}.sorted.tsv"
    echo "Starting bgzip"
    bgzip "~{atac_fragment_base}.sorted.tsv"
    echo "Starting tabix"
    tabix -s 1 -b 2 -e 3 -C "~{atac_fragment_base}.sorted.tsv.gz"
  >>>

  runtime {
    docker: docker_path
    disks: "local-disk ~{disk} HDD"
    memory: "${machine_mem_mb} MiB"
    cpu: nthreads
  }

  output {
    File gex_h5ad_file = "~{gex_base_name}.h5ad"
    File atac_h5ad_file = "~{atac_base_name}.h5ad"
    File atac_fragment_tsv = "~{atac_fragment_base}.sorted.tsv.gz"
    File atac_fragment_tsv_index = "~{atac_fragment_base}.sorted.tsv.gz.csi"
    File gex_whitelist_name_file = "gex_whitelist_used.txt"
    File atac_whitelist_name_file = "atac_whitelist_used.txt"
  }
}

task SlideseqH5adGeneration {

  input {
    #runtime values
    String warp_tools_docker_path
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
    Int machine_mem_mb = 32000
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
    docker: warp_tools_docker_path
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

task SingleNucleusSlideseqH5adOutput {

    input {
        #runtime values
        String warp_tools_docker_path
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
        description: "This task will converts some of the outputs of Optimus pipeline into a h5ad file"
    }

    parameter_meta {
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command <<<
        set -euo pipefail

        python3 /warptools/scripts/create_snrna_optimus_exons_h5ad.py \
        --annotation_file ~{annotation_file} \
        --cell_metrics ~{cell_metrics} \
        --gene_metrics ~{gene_metrics} \
        --count_matrix_1 ~{sparse_count_matrix} \
        --cell_id_1 ~{cell_id} \
        --gene_id_1  ~{gene_id} \
        --count_matrix_2 ~{sparse_count_matrix_exon} \
        --cell_id_2 ~{cell_id_exon} \
        --gene_id_2  ~{gene_id_exon} \
        --output_path_for_h5ad "~{input_id}" \
        --input_id ~{input_id} \
        ~{"--input_name " + input_name} \
        ~{"--input_id_metadata_field " + input_id_metadata_field} \
        ~{"--input_name_metadata_field " + input_name_metadata_field} \
        --expression_data_type "whole_transcript" \
        --pipeline_version ~{pipeline_version}
        
      >>>

    runtime {
        docker: warp_tools_docker_path
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

task SingleNucleusSmartSeq2H5adOutput {
    input {
        #runtime values
        String docker = "us.gcr.io/broad-gotc-prod/warp-tools:2.6.1"

        Array[File] alignment_summary_metrics
        Array[File] dedup_metrics
        Array[File] gc_bias_summary_metrics

        # introns counts
        Array[File] introns_counts
        # exons counts
        Array[File] exons_counts
        # annotation file
        File annotation_introns_added_gtf
        # name of the sample
        Array[String] input_ids
        Array[String]? input_names
        String? input_id_metadata_field
        String? input_name_metadata_field

        String pipeline_version
        Int preemptible = 3
        Int disk = 200
        Int machine_mem_mb = 8000
        Int cpu = 4
    }

    meta {
        description: "This task will convert output from the SmartSeq2SingleNucleus pipeline into a loom file. Contrary to the SmartSeq2 single cell where there is only RSEM counts, here we have intronic and exonic counts per gene name"
    }

    parameter_meta {
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command <<<
        set -euo pipefail

        declare -a introns_counts_files=(~{sep=' ' introns_counts})
        declare -a exons_counts_files=(~{sep=' ' exons_counts})
        declare -a output_prefix=(~{sep=' ' input_ids})
        declare -a alignment_summary_metrics_list=(~{sep=' 'alignment_summary_metrics})
        declare -a dedup_metrics_list=(~{sep=' 'dedup_metrics})
        declare -a gc_bias_summary_metrics_list=(~{sep=' 'gc_bias_summary_metrics})

        for (( i=0; i<${#introns_counts_files[@]}; ++i));
        do
        # creates a table with gene_id, gene_name, intron and exon counts
        echo "Running create_snss2_counts_csv."
        python /warptools/scripts/create_snss2_counts_csv.py \
        --in-gtf ~{annotation_introns_added_gtf} \
        --intron-counts ${introns_counts_files[$i]} \
        --exon-counts ${exons_counts_files[$i]}  \
        -o "${output_prefix[$i]}.exon_intron_counts.tsv"
        echo "Success create_snss2_counts_csv."

        # groups the QC file into one file
        echo "Running GroupQCs"
        GroupQCs -f "${alignment_summary_metrics_list[$i]}" "${dedup_metrics_list[$i]}" "${gc_bias_summary_metrics_list[$i]}" \
        -t Picard -o "${output_prefix[$i]}.Picard_group"
        echo "Success GroupQCs"

        # create the loom file
        echo "Running create_h5ad_snss2."
        python3 /warptools/scripts/create_h5ad_snss2.py \
        --qc_files "${output_prefix[$i]}.Picard_group.csv" \
        --count_results  "${output_prefix[$i]}.exon_intron_counts.tsv" \
        --output_h5ad_path "${output_prefix[$i]}" \
        --input_id ${output_prefix[$i]} \
        ~{"--input_id_metadata_field " + input_id_metadata_field} \
        ~{"--input_name_metadata_field " + input_name_metadata_field} \
        --pipeline_version ~{pipeline_version}

        echo "Success create_h5ad_snss2"
        done;
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{machine_mem_mb} MiB"
        disks: "local-disk ~{disk} HDD"
        disk: disk + " GB" # TES
        preemptible: preemptible
    }

    output {
        Array[File] h5ad_output = glob("*.h5ad")
        Array[File] exon_intron_counts = glob("*exon_intron_counts.tsv")
    }
}

task AggregateSmartSeq2H5ad {
    input {
        Array[File] h5ad_input
        String batch_id
        String pipeline_version
        String docker = "us.gcr.io/broad-gotc-prod/warp-tools:2.6.1"
        Int disk = 200
        Int machine_mem_mb = 4000
        Int cpu = 1
    }

    meta {
        description: "aggregate the H5AD output"
    }

    command {
        set -e

        # Merge the h5ad files
        python3 /warptools/scripts/ss2_h5ad_merge.py \
        --input-h5ad-files ~{sep=' ' h5ad_input} \
        --output-h5ad-file "~{batch_id}.h5ad" \
        --batch_id ~{batch_id} \
        --pipeline_version ~{pipeline_version}
    }

    output {
        File h5ad_output_file = "~{batch_id}.h5ad"
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
