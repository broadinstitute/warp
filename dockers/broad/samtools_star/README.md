# Samtools_STAR

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/samtools-star:1.0.0-1.11-2.7.10a-1642556627`
`

- __What is this image:__ This image is a lightweight alpine-based custom image for running Samtools, and the STAR aligner, it uses `us.gcr.io/broad-gotc-prod/samtools` as a base image.
- __What are Samtools, STAR:__  STAR is a fast universal RNA-seq aligner, [more info](https://github.com/alexdobin/STAR). See Samtools README for [more info](../samtools/README.md).
- __How to see tool version used in image:__ Please see below.

## Versioning

Samtools_STAR uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/samtools-star:1.0.0-1.11-2.7.10a-1642556627:<image-version>-<samtools-version>-<star-version>-<unix-timestamp>` 


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/samtools-star:1.0.0-1.11-2.7.10a-1642556627
$ docker inspect us.gcr.io/broad-gotc-prod/samtools-star:1.0.0-1.11-2.7.10a-1642556627
```

## Usage

### STAR

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/samtools-star:1.0.0-1.11-2.7.10a-1642556627 STAR
```

### Samtools

See Samtools README for [more info](../samtools/README.md)