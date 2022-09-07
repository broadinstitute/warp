# WARP RTools

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/rtools:1.0.0-1662498122`

- __What is this image:__ This image is a Debian-based custom image that contains R-based tools used in various WARP pipelines.
- __How to see tool version used in image:__ Please see below.

## Versioning

RTools uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/rtools:<image-version>-<unix-timestamp>`


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/rtools:1.0.0-1662498122
$ docker inspect us.gcr.io/broad-gotc-prod/rtools:1.0.0-1662498122
```

## Usage

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/pytools:1.0.0-1662498122
```