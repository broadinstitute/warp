# ATAC_BARCODES

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/atac_barcodes:1.0.4-1.0-1679935783`


- __What is this image:__ This image is a custom image for attaching cell barcodes from a 10x Read2 FASTQ to paired-read FASTQs.
- __What is ATAC_BARCODES:__ ATAC_BARCODES is a custom python script for attaching cell barcodes listed in one FASTQ to the readnames in two additional FASTQ files. 
- __How to see tool version used in image:__ Please see below.

## Versioning

ATAC_BARCODES image uses the following convention for verisoning:

#### `us.gcr.io/broad-gotc-prod/atac_barcodes:<image-version>-<atac_barcodes-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/atac_barcodes:1.0.4-1.0-1679935783
$ docker inspect us.gcr.io/broad-gotc-prod/atac_barcodes:1.0.4-1.0-1679935783
```

## Usage

### ATAC BARCODES 

```bash
$ docker run --rm -it \
    -v /files:/files \
    us.gcr.io/broad-gotc-prod/atac_barcodes:1.0.4-1.0-1679935783 bash
```