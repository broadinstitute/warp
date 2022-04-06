# RNA Seq fastp

## Quick reference

Copy and paste to pull this image

#### `us.gcr.io/broad-gotc-prod/fastp:1.0.0-0.20.1-1649253500`

- __What is this image:__ This image is a lightweight debian based image for running the fastp tool set within our RNA sequencing pipeline.
- __What is fastp:__ fastp from OpenGene is a tool designed to provide fast all-in-one preprocessing for FastQ files. See [here](https://github.com/OpenGene/fastp) for more information.
- __How to see tool version used in image:__ Please see below.

## Versioning

fastp uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/fastp:<image-version>-<fastp-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/fastp:1.0.0-0.20.1-1649253500
$ docker inspect us.gcr.io/broad-gotc-prod/fastp:1.0.0-0.20.1-1649253500
```

## Usage

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/fastp:1.0.0-0.20.1-1649253500 fastp
```