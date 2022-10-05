# Subread

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/subread:1.0.0-2.0.1-1662044537`


- __What is this image:__ This image is a lightweight alpine-based custom image for running the Subread sequencing data processing suite.
- __What is Subread:__  Subread contains a suite of high-performance software programs for processing next-generation sequencing data, [more info](http://subread.sourceforge.net).
- __How to see tool version used in image:__ Please see below.

## Versioning

Subread uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/subread:<image-version>-<subread-version>-<unix-timestamp>`


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/subread:1.0.0-2.0.1-1662044537
$ docker inspect us.gcr.io/broad-gotc-prod/subread:1.0.0-2.0.1-1662044537
```

## Usage

This image contains several tools. See [here](http://subread.sourceforge.net/subread.html) for a sample workflow.
To show the `featureCounts` help page, for example:

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/subread:1.0.0-2.0.1-1662044537 featureCounts
```