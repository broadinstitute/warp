version 1.0



task OptimusH5adGeneration {

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

task SingleNucleusOptimusH5adOutput {

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

    command {
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
        File h5ad_output = "~{input_id}.h5ad"
    }
}

task JoinMultiomeBarcodes {
    input {
    File atac_h5ad
    File atac_fragment
    File gex_h5ad
    File gex_whitelist
    File atac_whitelist

    Int nthreads = 1
    String cpuPlatform = "Intel Cascade Lake"
  }
    String gex_base_name = basename(gex_h5ad, ".h5ad")
    String atac_base_name = basename(atac_h5ad, ".h5ad")
    String atac_fragment_base = basename(atac_fragment, ".tsv")

    Int machine_mem_mb = ceil((size(atac_h5ad, "MiB") + size(gex_h5ad, "MiB") + size(atac_fragment, "MiB")) * 3) + 10000
    Int disk =  ceil((size(atac_h5ad, "GiB") + size(gex_h5ad, "GiB") + size(atac_fragment, "GiB")) * 5) + 10

  parameter_meta {
    atac_h5ad: "The resulting h5ad from the ATAC workflow."
    atac_fragment: "The resulting fragment TSV from the ATAC workflow."
    gex_h5ad: "The resulting h5ad from the Optimus workflow."
    gex_whitelist: "Whitelist used for gene expression barcodes."
    atac_whitelist: "Whitelist used for ATAC barcodes."
  }

  command <<<
    set -e pipefail

    python3 <<CODE

    # set parameters
    atac_h5ad = "~{atac_h5ad}"
    atac_fragment = "~{atac_fragment}"
    gex_h5ad = "~{gex_h5ad}"
    gex_whitelist = "~{gex_whitelist}"
    atac_whitelist = "~{atac_whitelist}"

    # import anndata to manipulate h5ad files
    import anndata as ad
    import pandas as pd
    print("Reading ATAC h5ad:")
    print("~{atac_h5ad}")
    print("Read ATAC fragment file:")
    print("~{atac_fragment}")
    print("Reading Optimus h5ad:")
    print("~{gex_h5ad}")
    atac_data = ad.read_h5ad("~{atac_h5ad}")
    gex_data = ad.read_h5ad("~{gex_h5ad}")
    atac_tsv = pd.read_csv("~{atac_fragment}", sep="\t", names=['chr','start', 'stop', 'barcode','n_reads'])
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
    # write out the files
    gex_data.write("~{gex_base_name}.h5ad")
    atac_data.write_h5ad("~{atac_base_name}.h5ad")
    df_fragment.to_csv("~{atac_fragment_base}.tsv", sep='\t', index=False, header = False)
    CODE
    # sorting the file
    echo "Sorting file"
    sort -k1,1V -k2,2n "~{atac_fragment_base}.tsv" > "~{atac_fragment_base}.sorted.tsv"
    echo "Starting bgzip"
    bgzip "~{atac_fragment_base}.sorted.tsv"
    echo "Starting tabix"
    tabix -s 1 -b 2 -e 3 "~{atac_fragment_base}.sorted.tsv.gz"

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/snapatac2:1.0.4-2.3.1-1700590229"
    disks: "local-disk ~{disk} HDD"
    memory: "${machine_mem_mb} MiB"
    cpu: nthreads
  }

  output {
    File gex_h5ad_file = "~{gex_base_name}.h5ad"
    File atac_h5ad_file = "~{atac_base_name}.h5ad"
    File atac_fragment_tsv = "~{atac_fragment_base}.sorted.tsv.gz"
    File atac_fragment_tsv_tbi = "~{atac_fragment_base}.sorted.tsv.gz.tbi"
  }
}
