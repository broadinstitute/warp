# Samtools_BWA

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091`


- __What is this image:__ This image is a lightweight alpine-based custom image for running Samtools and BWA, it uses `us.gcr.io/broad-gotc-prod/samtools` as a base image.
- __What are Samtools and BWA:__ BWA is a software package for mapping DNA sequewnces against a large reference genome, [more info](https://github.com/lh3/bwa). Samtools is a suite of programs for interacting with high-throughput sequencing data. See [here](https://github.com/samtools/samtools) more information. 
- __How to see tool version used in image:__ Please see below.

## Versioning

Samtools_BWA uses the following convention for verisoning:

#### `us.gcr.io/broad-gotc-prod/samtools-bwa:<image-version>-<bwa-version>-<unix-timestamp>` 

To see the Samtools version being used please reference the base image [here](..samtools/README.md).

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091
$ docker inspect us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091
```

## Usage

### BWA 

```bash
$ docker run --rm -it \
    -v /bwa-files:bwa-files \
    us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091 /usr/gitc/bwa mem \
    /bwa-files/ref.fa /bwa-files/reads.fq > /bwa-files/aln.sam
```

### Samtools

See Samtools README for [more info](../samtools/README.md)
