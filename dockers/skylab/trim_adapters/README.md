# TrimAdapters

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/trim-adapters:1.0.0-1.04.807-1659980243`
`

- __What is this image:__ This image is an Ubuntu-based custom image for running TrimAdapters.
- __What is TrimAdapters:__  TrimAdapters includes tools for processing biological sequencing data, [more info](http://expressionanalysis.github.io/ea-utils/).
- __How to see tool version used in image:__ Please see below.

## Versioning

TrimAdapters uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/trim-adapters:<image-version>-<star-version>-<unix-timestamp>`


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/trim-adapters:1.0.0-1.04.807-1659980243
$ docker inspect us.gcr.io/broad-gotc-prod/trim-adapters:1.0.0-1.04.807-1659980243
```

## Usage

Show the `fastq-mcf` help screen:

```bash
$ docker run --rm -it \
    us.gcr.io/trim-adapters:1.0.0-1.04.807-1659980243 fastq-mcf -h
```