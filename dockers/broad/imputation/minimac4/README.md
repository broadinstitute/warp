# Imputation Minimac4

## Quick reference

Copy and paste to pull this image

#### `us.gcr.io/broad-gotc-prod/imputation-minimac4:1.0.6-1.0.2-1663948783`

- __What is this image:__ This image is a lightweight alpine-based image for running Minimac4 in the [Imputation pipeline](../../../../pipelines/broad/arrays/imputation/Imputation.wdl).
- __What is Minimac4:__ Minimac4 is a low-memory and computationally efficient piece of software for genotype imputation. See [here](https://github.com/statgen/Minimac4) more information.
- __How to see Minimac4 version used in image:__ Please see below.

## Versioning

The Imputation Minimac4 image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/samtools:<image-version>-<minimac4-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/imputation-minimac4:1.0.6-1.0.2-1663948783
$ docker inspect us.gcr.io/broad-gotc-prod/imputation-minimac4:1.0.6-1.0.2-1663948783
```

## Usage

### Display default menu

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/imputation-minimac4:1.0.6-1.0.2-1663948783 /usr/gitc/minimac4
```