# WARP Python Tools

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/pytools:1.0.0-1661263730`

- __What is this image:__ This image is a Debian-based custom image that contains Python scripts used in various WARP pipelines.
- __How to see tool version used in image:__ Please see below.

## Versioning

PyTools uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/pytools:<image-version>-<unix-timestamp>`


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/pytools:1.0.0-1661263730
$ docker inspect us.gcr.io/broad-gotc-prod/pytools:1.0.0-1661263730
```

## Usage

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/pytools:1.0.0-1661263730 <script.py>
```

## Scripts

This image contains the following scripts:

* `breakoutSnap.py` extracts the data in a snap file as csv files
* `create-merged-npz-output.py` takes a barcode.tsv, feature.tsv and matrix.mtx from STAR alignment outputs and creates 2 npy files and an npz file for row_index, col_index and the matrix. These files are required in the empty_drop step.
* `create_snss2_counts_csv.py` creates a csv file containing intron and exon counts from the Single Nucleus Smart-Seq2 pipeline
* `loomCompare.py` compares differences between loom files
* `ss2_loom_merge.py` creates a single loom file from multiple single sample loom files
* `makeCompliantBAM.py` make a BAM file with cellular barcodes in the read names compliant by moving them to the CB tag

The following scripts create a loom file from counts, metadata, and metrics from each pipeline:
* `create_loom_optimus.py` for Optimus pipeline
* `create_loom_snss2.py` for Single Nucleus Smart-Seq2 pipeline
* `create_snrna_optimus.py` for Optimus in `sn_rna` mode with `count_exons=false`
* `create_snrna_optimus_counts.py` for Optimus in `sn_rna` mode with `count_exons=true`

