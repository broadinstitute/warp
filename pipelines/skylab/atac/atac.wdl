version 1.0

# This version of ATAC is creating a fragment file from the aligned bam file. This is NOT intended to be merged into WARP.

workflow ATAC {
  meta {
    description: "Processing for single-cell ATAC-seq data from the level of raw fastq reads. This is the first step of the multiome pipeline. ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a technique used in molecular biology to assess genome-wide chromatin accessibility. This pipeline processes 10x Genomics Multiome ATAC FASTQ files."
    allowNestedInputs: true
  }

  input {
    # Bam file to create fragment file from
    File aligned_bam
    # Text file containing chrom_sizes for genome build (i.e. hg38)
    File chrom_sizes
    #File for annotations for calculating ATAC TSSE
    File annotations_gtf
    # Option for running files with preindex
    Boolean preindex = false
    # Additional library aliquot ID
    String? atac_nhash_id
    #Expected cells from library preparation
    Int atac_expected_cells = 3000
    # Output prefix/base name for all intermediate files and pipeline outputs
    String input_id
  }

  String pipeline_version = "n/a"

  call CreateFragmentFile {
    input:
      bam = aligned_bam,
      chrom_sizes = chrom_sizes,
      annotations_gtf = annotations_gtf,
      preindex = preindex,
      docker_path = "us.gcr.io/broad-gotc-prod/snapatac2:2.0.0",
      atac_nhash_id = atac_nhash_id,
      atac_expected_cells = atac_expected_cells,
      input_id = input_id
  }

  output {
    File bam_aligned_output = aligned_bam
    File fragment_file = CreateFragmentFile.fragment_file
    File fragment_file_index = CreateFragmentFile.fragment_file_index
    File snap_metrics = CreateFragmentFile.Snap_metrics
    File library_metrics_file = CreateFragmentFile.atac_library_metrics
  }
}

# This task creates a fragment file from the aligned bam file. It also calculates the TSSE metrics and library metrics for the fragment file.
task CreateFragmentFile {
  input {
    File bam
    File annotations_gtf
    File chrom_sizes
    Boolean preindex
    Int disk_size = 500
    Int mem_size = 64
    Int nthreads = 4
    String cpuPlatform = "Intel Cascade Lake"
    String docker_path
    String atac_nhash_id = ""
    String input_id
    Int atac_expected_cells = 3000
    String gtf_path = annotations_gtf
    Array[String] mito_list = ['chrM', 'M']
  }

  parameter_meta {
    bam: "Aligned bam with CB in CB tag. This is the output of the BWAPairedEndAlignment task."
    chrom_sizes: "Text file containing chrom_sizes for genome build (i.e. hg38)."
    annotations_gtf: "GTF for SnapATAC2 to calculate TSS sites of fragment file."
    disk_size: "Disk size used in create fragment file step."
    mem_size: "The size of memory used in create fragment file."
    docker_path: "The docker image path containing the runtime environment for this task"
  }

  command <<<
    set -euo pipefail
    set -x

    python3 <<CODE

    import snapatac2.preprocessing as pp
    import snapatac2 as snap
    import scanpy as sc
    import numpy as np
    import polars as pl
    import anndata as ad
    from collections import OrderedDict
    import csv

    # set parameters
    bam = "~{bam}"
    input_id = "~{input_id}"
    chrom_sizes = "~{chrom_sizes}"
    atac_gtf = "~{annotations_gtf}"
    preindex = "~{preindex}"
    atac_nhash_id = "~{atac_nhash_id}"
    expected_cells = ~{atac_expected_cells}
    mito_list = "~{sep=' ' mito_list}"

    print(mito_list)
    mito_list = mito_list.split(" ")
    print("Mitochondrial chromosomes:", mito_list)

    # calculate chrom size dictionary based on text file
    chrom_size_dict={}
    with open('~{chrom_sizes}', 'r') as f:
      for line in f:
        key, value = line.strip().split()
        chrom_size_dict[str(key)] = int(value)

    # use snap atac2
    import snapatac2.preprocessing as pp
    import snapatac2 as snap
    import scanpy as sc
    import numpy as np
    import polars as pl
    import anndata as ad
    from collections import OrderedDict
    import csv

    # extract CB or BB (if preindex is true) tag from bam file to create fragment file
    if preindex == "true":
      data = pp.recipe_10x_metrics("~{bam}", "~{input_id}.fragments.tsv", "temp_metrics.h5ad", is_paired=True, barcode_tag="BB", chrom_sizes=chrom_size_dict, gene_anno=atac_gtf, peaks=None, chrM=mito_list)
    elif preindex == "false":
      data = pp.recipe_10x_metrics("~{bam}", "~{input_id}.fragments.tsv", "temp_metrics.h5ad", is_paired=True, barcode_tag="CB", chrom_sizes=chrom_size_dict, gene_anno=atac_gtf, peaks=None, chrM=mito_list)

    # Add NHashID to metrics
    data = OrderedDict({'NHashID': atac_nhash_id, **data})

    # Calculate atac percent target
    print("Calculating percent target")
    number_of_cells = data['Cells']['Number_of_cells']
    print("Print number of cells", number_of_cells)
    atac_percent_target = number_of_cells / expected_cells*100
    print("Setting percent target in nested dictionary")
    data['Cells']['atac_percent_target'] = atac_percent_target

    # Flatten the dictionary
    flattened_data = []
    for category, metrics in data.items():
        if isinstance(metrics, dict):
            for metric, value in metrics.items():
                flattened_data.append((metric, value))
        else:
            flattened_data.append((category, metrics))

    # Convert the flattened keys to lowercase (except for 'NHashID')
    flattened_data = [(metric if metric == 'NHashID' else str(metric).lower(), value) for metric, value in flattened_data]

    # Write to CSV
    csv_file_path = "~{input_id}_~{atac_nhash_id}_library_metrics.csv"
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(flattened_data)  # Write data

    print(f"Dictionary successfully written to {csv_file_path}")

    atac_data = ad.read_h5ad("temp_metrics.h5ad")
    # Add nhash_id to h5ad file as unstructured metadata
    atac_data.uns['NHashID'] = atac_nhash_id

    # Add GTF to uns field
    # Original path from args.annotation_file
    gtf_path = "~{gtf_path}"  # e.g., 'gs://gcp-public-data--broad-references/hg38/v0/star/v2_7_10a/modified_v43.annotation.gtf'

    atac_data.uns["reference_gtf_file"] = gtf_path
    # calculate tsse metrics
    snap.metrics.tsse(atac_data, atac_gtf)
    # Write new atac file
    atac_data.write_h5ad("~{input_id}.metrics.h5ad")

    CODE

    # sorting the file
    echo "Sorting file"
    sort -k1,1V -k2,2n "~{input_id}.fragments.tsv" > "~{input_id}.fragments.sorted.tsv"
    echo "Starting bgzip"
    bgzip "~{input_id}.fragments.sorted.tsv"
    echo "Starting tabix"
    tabix -s 1 -b 2 -e 3 -C "~{input_id}.fragments.sorted.tsv.gz"
  >>>

  runtime {
    docker: docker_path
    disks: "local-disk ${disk_size} SSD"
    memory: "${mem_size} GiB"
    cpu: nthreads
    cpuPlatform: cpuPlatform
  }

  output {
    File fragment_file = "~{input_id}.fragments.sorted.tsv.gz"
    File fragment_file_index = "~{input_id}.fragments.sorted.tsv.gz.csi"
    File Snap_metrics = "~{input_id}.metrics.h5ad"
    File atac_library_metrics = "~{input_id}_~{atac_nhash_id}_library_metrics.csv"
  }
}