# BWA

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/bwa:1.0.0-0.7.17-1660770463`


- __What is this image:__ This image is an Ubuntu-based custom image for running BWA.
- __What is BWA:__ BWA is a software package for mapping DNA sequences against a large reference genome, [more info](http://bio-bwa.sourceforge.net).
- __How to see tool version used in image:__ Please see below.

## Versioning

BWA uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/bwa:<image-version>-<bwa-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/bwa:1.0.0-0.7.17-1660770463
$ docker inspect us.gcr.io/broad-gotc-prod/bwa:1.0.0-0.7.17-1660770463
```

## Usage

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/bwa:1.0.0-0.7.17-1660770463 \
    bwa
```