# Snaptools_BWA

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/snaptools-bwa:1.0.0-1.4.8-0.7.17-1660844602`

- __What is this image:__ This image is a Debian-based custom image for running Snaptools and BWA.
- __What are Snaptools and BWA:__ BWA is a software package for mapping DNA sequences against a large reference genome, [more info](http://bio-bwa.sourceforge.net). Snaptools is a Python module for working with snap files, [more info](https://github.com/r3fang/SnapTools).
- __How to see tool version used in image:__ Please see below.

## Versioning

Snaptools_BWA uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/snaptools_bwa:<image-version>-<snaptools_version>-<bwa-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/snaptools-bwa:1.0.0-1.4.8-0.7.17-1660844602
$ docker inspect us.gcr.io/broad-gotc-prod/snaptools-bwa:1.0.0-1.4.8-0.7.17-1660844602
```

## Usage

### BWA 

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/snaptools-bwa:1.0.0-1.4.8-0.7.17-1660844602 \
    bwa
```

### Snaptools

See Snaptools GitHub for [more info](https://github.com/r3fang/SnapTools).