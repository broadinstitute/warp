# STAR

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658334187`

- __What is this image:__ This image is a lightweight alpine-based custom image for running the STAR aligner.
- __What is STAR:__  STAR is a fast universal RNA-seq aligner, [more info](https://github.com/alexdobin/STAR).
- __How to see tool version used in image:__ Please see below.

## Versioning

STAR uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/star:<image-version>-<star-version>-<unix-timestamp>`


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658334187
$ docker inspect us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658334187
```

## Usage

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658334187 STAR
```