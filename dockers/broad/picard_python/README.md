# Picard_Python

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/picard-python:1.0.0-2.26.10-1663951039`

- __What is this image:__ This image is a lightweight alpine-based custom image for running Picard and Python, it uses `python:3.8-alpine` as a base image.
- __What is Picard:__ Picard is a set of command line tools for manipulating high-throughput sequencing (HTS) data and formats, [more info](https://github.com/broadinstitute/picard).
- __How to see tool version used in image:__ Please see below.

## Versioning

Picard_Python uses the following convention for verisoning:

#### `us.gcr.io/broad-gotc-prod/picard-python:1.0.0-2.26.10-1663951039` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/picard-python:1.0.0-2.26.10-1663951039     
$ docker inspect us.gcr.io/broad-gotc-prod/picard-python:1.0.0-2.26.10-1663951039 
```

## Usage

### Picard

See Picard GitHub for [more info](https://github.com/broadinstitute/picard)
