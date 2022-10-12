# HISAT2

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/hisat2:1.0.0-1662998171`

- __What is this image:__ This image is an Ubuntu-based custom image that contains the HISAT2 alignment program.
- __What is HISAT2:__ HISAT2 is a fast and sensitive spliced alignment program for mapping RNA-seq reads, [more info](https://ccb.jhu.edu/software/hisat/manual.shtml#running-hisat).
- __How to see tool version used in image:__ Please see below.

## Versioning

HISAT2 uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/hisat2:<image-version>-<unix-timestamp>`


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/hisat2:1.0.0-1662998171
$ docker inspect us.gcr.io/broad-gotc-prod/hisat2:1.0.0-1662998171
```

## Usage

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/hisat2:1.0.0-1662998171 hisat2
```