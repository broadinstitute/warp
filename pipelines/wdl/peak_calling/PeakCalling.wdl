version 1.0

import "../../../tasks/wdl/Utilities.wdl" as utils

workflow PeakCalling {

  meta {
    allowNestedInputs: true
  }

  input {

      File annotations_gtf
      File metrics_h5ad
      File chrom_sizes
      String output_base_name

      # SnapATAC2 parameters
      Int min_counts = 5000
      Int min_tsse = 10
      Int max_counts = 100000
      Float probability_threshold = 0.5

      # Runtime attributes/docker
      Int disk_size = 500
      Int mem_size = 64
      Int nthreads = 4

      String cloud_provider
  }

  String pipeline_version = "1.0.1"

  # Determine docker prefix based on cloud provider
  String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
  String acr_docker_prefix = "dsppipelinedev.azurecr.io/"
  String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

  # Docker image names
  String snap_atac_docker = "snapatac2:2.0.0"

  # Make sure either 'gcp' or 'azure' is supplied as cloud_provider input. If not, raise an error
  if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
      call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
        input:
        message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
    }
  }

  call PeakCalling {
      input:
        annotations_gtf = annotations_gtf,
        metrics_h5ad = metrics_h5ad,
        chrom_sizes = chrom_sizes,
        output_base_name = output_base_name,
        docker_path = docker_prefix + snap_atac_docker

  }

    output {
        File cellbybin_h5ad = PeakCalling.cellbybin_h5ad
        File cellbypeak_h5ad = PeakCalling.cellbypeak_h5ad
        String pipeline_version_out = pipeline_version
    }
}

# peak calling using SnapATAC2
task PeakCalling {
    input {
        File annotations_gtf
        File metrics_h5ad
        File chrom_sizes
        String output_base_name

        # SnapATAC2 parameters
        Int min_counts = 5000
        Int min_tsse = 10
        Int max_counts = 100000
        Float probability_threshold = 0.5

        # Runtime attributes/docker
        String docker_path
        Int disk_size = 500
        Int mem_size = 64
        Int nthreads = 4
    }

    parameter_meta {
        annotations_gtf: "GTF for SnapATAC2 to calculate TSS sites of fragment file."
        disk_size: "Disk size used in create fragment file step."
        mem_size: "The size of memory used in create fragment file."
        docker_path: "The docker image path containing the runtime environment for this task"
    }

    command <<<
        set -euo pipefail
        set -x

        python3 <<CODE

        # use snap atac2
        import snapatac2 as snap
        import scanpy as sc
        import numpy as np
        import polars as pl
        import pandas as pd

        output_base_name = "~{output_base_name}"
        atac_gtf = "~{annotations_gtf}"
        metrics_h5ad = "~{metrics_h5ad}"
        chrom_sizes = "~{chrom_sizes}"
        min_counts = "~{min_counts}"
        min_tsse = "~{min_tsse}"
        max_counts = "~{max_counts}"
        probability_threshold = "~{probability_threshold}"

        probability_threshold = float(probability_threshold)

        print("Peak calling starting...")
        atac_data = snap.read(metrics_h5ad)

        # Calculate and plot the size distribution of fragments
        print("Calculating fragment size distribution")
        snap.pl.frag_size_distr(atac_data)
        print(atac_data)

        # Filter cells
        print("Filtering cells")
        snap.pp.filter_cells(atac_data, min_counts=min_counts, min_tsse=min_tsse, max_counts=max_counts)
        print(atac_data)

        # Create a cell by bin matrix containing insertion counts across genome-wide 500-bp bins.
        print("Creating cell by bin matrix")
        atac_data_mod = snap.pp.add_tile_matrix(atac_data, inplace=False)
        print("set obsm")
        atac_data_mod.obsm["fragment_paired"] =  atac_data.obsm["fragment_paired"]
        print("set all uns")
        for key in atac_data.uns.keys():
          print("set ",key)
          atac_data_mod.uns[key] = atac_data.uns[key]
        print(atac_data_mod)

        # Feature selection
        print("Feature selection")
        snap.pp.select_features(atac_data_mod)
        print(atac_data_mod)

        # Run customized scrublet algorithm to identify potential doublets
        print("Run scrublet to identify potential doublets")
        snap.pp.scrublet(atac_data_mod)
        print(atac_data_mod)

        # Employ spectral embedding for dimensionality reduction
        print("Employ spectral embedding for dimensionality reduction")
        snap.tl.spectral(atac_data_mod)
        print(atac_data_mod)

        # Filter doublets based on scrublet scores
        print("Filter doublets based on scrublet scores")
        snap.pp.filter_doublets(atac_data_mod, probability_threshold=probability_threshold)
        print(atac_data_mod)

        # Check if the matrix is empty
        if atac_data_mod.n_obs == 0:
          raise ValueError("Matrix is empty after filtering doublets: Try increasing the probability_threshold.")

        # Perform graph-based clustering to identify cell clusters.
        # Build a k-nearest neighbour graph using snap.pp.knn
        print("Perform knn graph-based clustering to identify cell clusters")
        snap.pp.knn(atac_data_mod)
        print(atac_data_mod)

        # Use the Leiden community detection algorithm to identify densely-connected subgraphs/clusters in the graph
        print("Use the Leiden community detection algorithm to identify densely-connected subgraphs/clusters in the graph")
        snap.tl.leiden(atac_data_mod)
        print(atac_data_mod)

        # Create the cell by gene activity matrix
        print("Create the cell by gene activity matrix")
        gene_mat = snap.pp.make_gene_matrix(atac_data_mod, gene_anno=atac_gtf)
        print(atac_data_mod)

        # Normalize the gene matrix
        print("Normalize the gene matrix")
        gene_mat.obs['leiden'] = atac_data_mod.obs['leiden']
        sc.pp.normalize_total(gene_mat)
        sc.pp.log1p(gene_mat)
        sc.tl.rank_genes_groups(gene_mat, groupby="leiden", method="wilcoxon")

        for i in np.unique(gene_mat.obs['leiden']):
            markers = sc.get.rank_genes_groups_df(gene_mat, group=i).head(7)['names']
            print(f"Cluster {i}: {', '.join(markers)}")

        print("Peak calling using MACS3")
        snap.tl.macs3(atac_data_mod, groupby='leiden', n_jobs=1)

        print("Merge peaks and create peak matrix")
        # read chrom sizes
        chromsize_dict = pd.read_csv(chrom_sizes, sep='\t', header=None)
        chromsize_dict = pd.Series(chromsize_dict[1].values, index=chromsize_dict[0]).to_dict()
        # merge peaks and create peak matrix
        peaks = snap.tl.merge_peaks(atac_data_mod.uns['macs3'], chromsize_dict)
        peak_matrix = snap.pp.make_peak_matrix(atac_data_mod, use_rep=peaks['Peaks'])

        print("Convert pl.DataFrame to pandas DataFrame")
        # Convert pl.DataFrame to pandas DataFrame
        for key in atac_data_mod.uns.keys():
          if isinstance(atac_data_mod.uns[key], pl.DataFrame):
              print(key)
              atac_data_mod.uns[key] = atac_data_mod.uns[key].to_pandas()

        print("Write into h5ad file")
        atac_data_mod.write_h5ad("~{output_base_name}.cellbybin.h5ad")
        peak_matrix.write_h5ad("~{output_base_name}.cellbypeak.h5ad")

        CODE
    >>>

    runtime {
        docker: docker_path
        disks: "local-disk ${disk_size} SSD"
        memory: "${mem_size} GiB"
        cpu: nthreads
    }

    output {
        File cellbybin_h5ad = "~{output_base_name}.cellbybin.h5ad"
        File cellbypeak_h5ad = "~{output_base_name}.cellbypeak.h5ad"
    }
}
