# Imputation Eagle

## Quick reference

Copy and paste to pull this image

#### `us.gcr.io/broad-gotc-prod/imputation-eagle:1.0.0-2.4-1633695564`

- __What is this image:__ This image is a debian-based image for running Eagle in the [Imputation pipeline](../../../../pipelines/broad/arrays/imputation/Imputation.wdl).
- __What is Eagle:__ Eagle is an algorithm for estimating haplotype phases either wihin a genotyped cohort or using a phased reference panel. See [here](https://github.com/poruloh/Eagle) more information.
- __How to see Eagle version used in image:__ Please see below.

## Versioning

The Imputation Eagle image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/imputation-eagle:<image-version>-<eagle-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/imputation-eagle:1.0.0-2.4-1633695564
$ docker inspect us.gcr.io/broad-gotc-prod/imputation-eagle:1.0.0-2.4-1633695564
```

## Usage

### Display default menu

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/imputation-eagle:1.0.0-2.4-1633695564 /usr/gitc/eagle 
```