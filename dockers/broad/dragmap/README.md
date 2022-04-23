# Dragmap

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/dragmap:1.1.2-1.2.1-2.26.10-1.11-1643839530`

- __What is this image:__ This image is a RHEL-based custom image for running DRAGMAP, Picard and SAMTOOLS, it uses `centos:8` as a base image.
- __What are Dragmap, Picard and Samtools:__ Dragmap is the Dragen mapper/aligner Open Source Software, [more info](https://github.com/Illumina/DRAGMAP). Picard is a set of command line tools for manipulating high-throughput sequencing (HTS) data and formats, [more info](https://github.com/broadinstitute/picard). Samtools is a suite of programs for interacting with high-throughput sequencing data. See [here](https://github.com/samtools/samtools) more information.
- __How to see tool version used in image:__ Please see below.

## Versioning

Dragmap uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/dragmap:<image-version>-<dragmap-version>-<picard-version>-<samtools-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/dragmap:1.1.2-1.2.1-2.26.10-1.11-1643839530
$ docker inspect us.gcr.io/broad-gotc-prod/dragmap:1.1.2-1.2.1-2.26.10-1.11-1643839530
```

## Usage

### Dragmap 

See Dragmap GitHub for [more info](https://github.com/Illumina/DRAGMAP)

### Picard

See Picard GitHub for [more info](https://github.com/broadinstitute/picard)

### Samtools

```bash
$ docker run --rm -it \
    -v /bamfiles:/bamfiles \
    us.gcr.io/broad-gotc-prod/dragmap:1.1.2-1.2.1-2.26.10-1.11-1643839530 samtools view -H /bamfiles/<bam-file>
```