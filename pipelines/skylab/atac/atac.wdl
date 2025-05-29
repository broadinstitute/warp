version 1.0

workflow ATAC {
  input {
    File fragment_file
    File chrom_sizes
    File annotations_gtf
    String? atac_nhash_id
    Int atac_expected_cells = 3000
    String input_id
  }

  String pipeline_version = "n/a"

  # Call GenerateAtacMetrics with the correct fragment file
  call GenerateAtacMetrics {
    input:
      fragment_file = fragment_file,
      chrom_sizes = chrom_sizes,
      annotations_gtf = annotations_gtf,
      docker_path = "us.gcr.io/broad-gotc-prod/snapatac2:2.0.0",
      atac_nhash_id = atac_nhash_id,
      atac_expected_cells = atac_expected_cells,
      input_id = input_id
  }

  output {
    File Snap_metrics = GenerateAtacMetrics.Snap_metrics
    File library_metrics_file = GenerateAtacMetrics.atac_library_metrics
  }
}

  task GenerateAtacMetrics {
  input {
    File fragment_file
    File annotations_gtf
    File chrom_sizes
    Int disk_size = 500
    Int mem_size = 64
    Int nthreads = 4
    String cpuPlatform = "Intel Cascade Lake"
    String docker_path
    String atac_nhash_id = ""
    String input_id
    Int atac_expected_cells = 3000
  }

  command <<<
    set -euo pipefail
    set -x

python3 <<CODE
import snapatac2 as snap
import scanpy as sc
import anndata as ad
import csv
from collections import OrderedDict

# Import fragment file into AnnData
adata = snap.pp.import_data(
    fragment_file="~{fragment_file}",
    chrom_sizes="~{chrom_sizes}",
    is_paired=True,
    barcode_tag="CB"  # TODO or "BB" if your barcodes are indexed
)

# Calculate TSSE
snap.metrics.tsse(adata, "~{annotations_gtf}")

# Add metadata
adata.uns["NHashID"] = "~{atac_nhash_id}"
adata.uns["reference_gtf_file"] = "~{annotations_gtf}"

# Write H5AD
adata.write_h5ad("~{input_id}.metrics.h5ad")

# Generate metrics
number_of_cells = adata.n_obs
atac_percent_target = number_of_cells / ~{atac_expected_cells} * 100

# Write metrics CSV
metrics = OrderedDict({
    "NHashID": "~{atac_nhash_id}",
    "number_of_cells": number_of_cells,
    "atac_percent_target": atac_percent_target
})

with open("~{input_id}_~{atac_nhash_id}_library_metrics.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(metrics.items())
CODE
  >>>

  runtime {
    docker: docker_path
    disks: "local-disk ${disk_size} SSD"
    memory: "${mem_size} GiB"
    cpu: nthreads
    cpuPlatform: cpuPlatform
  }

  output {
    File Snap_metrics = "~{input_id}.metrics.h5ad"
    File atac_library_metrics = "~{input_id}_~{atac_nhash_id}_library_metrics.csv"
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

    # set parameters
    bam = "~{bam}"
    input_id = "~{input_id}"
    chrom_sizes = "~{chrom_sizes}"
    atac_gtf = "~{annotations_gtf}"
    preindex = "~{preindex}"
    atac_nhash_id = "~{atac_nhash_id}"
    expected_cells = ~{atac_expected_cells}

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
      data = pp.recipe_10x_metrics("~{bam}", "~{input_id}.fragments.tsv", "temp_metrics.h5ad", is_paired=True, barcode_tag="BB", chrom_sizes=chrom_size_dict, gene_anno=atac_gtf, peaks=None)
    elif preindex == "false":
      data = pp.recipe_10x_metrics("~{bam}", "~{input_id}.fragments.tsv", "temp_metrics.h5ad", is_paired=True, barcode_tag="CB", chrom_sizes=chrom_size_dict, gene_anno=atac_gtf, peaks=None)

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