# Expression Analysis Utilities (EA-Utils)

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/ea-utils:1.0.0-1.04.807-1659990665`

- __What is this image:__ This image is an Ubuntu-based custom image for running ea-utils.
- __What is ea-utils:__  ea-utils includes tools for processing biological sequencing data, [more info](http://expressionanalysis.github.io/ea-utils/).
- __How to see tool version used in image:__ Please see below.

## Versioning

ea-utils uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/ea-utils:<image-version>-<star-version>-<unix-timestamp>`


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/ea-utils:1.0.0-1.04.807-1659990665
$ docker inspect us.gcr.io/broad-gotc-prod/ea-utils:1.0.0-1.04.807-1659990665
```

## Usage

Show the `fastq-mcf` help screen:

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/ea-utils:1.0.0-1.04.807-1659990665 fastq-mcf -h
```