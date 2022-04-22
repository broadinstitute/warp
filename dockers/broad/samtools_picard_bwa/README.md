# Samtools_Picard_BWA

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748`
`

- __What is this image:__ This image is a lightweight alpine-based custom image for running Samtools, BWA and Picard, it uses `us.gcr.io/broad-gotc-prod/samtools` as a base image.
- __What are Samtools, BWA and Picard:__ BWA is a software package for mapping DNA sequewnces against a large reference genome, [more info](https://github.com/lh3/bwa). Picard is a set of command line tools for manipulating high-throughput sequencing (HTS) data and formats, [more info](https://github.com/broadinstitute/picard). See Samtools README for [more info](../samtools/README.md).
- __How to see tool version used in image:__ Please see below.

## Versioning

Samtools_Picard_BWA uses the following convention for verisoning:

#### `us.gcr.io/broad-gotc-prod/samtools_picard_bwa:<image-version>-<bwa-version>-<picard-version>-<unix-timestamp>` 

To see the Samtools version being used please reference the base image [here](..samtools/README.md)

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748
$ docker inspect us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748
```

## Usage

### BWA 

```bash
$ docker run --rm -it \
    -v /bwa-files:bwa-files \
    us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748 /usr/gitc/bwa mem \
    /bwa-files/ref.fa /bwa-files/reads.fq > /bwa-files/aln.sam
```

### Picard

See Picard GitHub for [more info](https://github.com/broadinstitute/picard)

### Samtools

See Samtools README for [more info](../samtools/README.md)