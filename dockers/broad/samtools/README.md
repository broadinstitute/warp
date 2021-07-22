# Samtools

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616`

- __What is this image:__ This image is a lightweight alpine-based image for running SAMTOOLS.
- __What is Samtools:__ Samtools is a suite of programs for interacting with high-throughput sequencing data. See [here](https://github.com/samtools/samtools) more information.
- __How to see Samtools version used in image:__ Please see below.

## Versioning

The Samtools image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/samtools:<image-version>-<samtools-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616
$ docker inspect us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616
```

## Usage

### Display default menu

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616 samtools 
```

### Viewing a BAM file

```bash
$ docker run --rm -it \
    -v /bamfiles:/bamfiles \
    us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616 samtools view -H /bamfiles/<bam-file>
```