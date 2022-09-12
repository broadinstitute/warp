# RSEM

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/rsem:1.0.0-1663016024`

- __What is this image:__ This image is an Ubuntu-based custom image that contains the RSEM tool suite.
- __What is RSEM:__ RSEM is a software package for estimating gene and isoform expression levels from RNA-Seq data, [more info](http://deweylab.biostat.wisc.edu/rsem/).
- __How to see tool version used in image:__ Please see below.

## Versioning

RSEM uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/rsem:<image-version>-<unix-timestamp>`


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/rsem:1.0.0-1663016024
$ docker inspect us.gcr.io/broad-gotc-prod/rsem:1.0.0-1663016024
```

## Usage

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/rsem:1.0.0-1663016024 rsem-prepare-reference --help
```