# SNAPATAC2

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/snapatac2:1.0.0-2.2.0-1679070508`


- __What is this image:__ This image is a custom image for running SnapATAC2.
- __What SanpATAC2:__ SnatATAC2 is a software package for , [more info](). 
- __How to see tool version used in image:__ Please see below.

## Versioning

SnapATAC2 image uses the following convention for verisoning:

#### `us.gcr.io/broad-gotc-prod/snapatac2:<image-version>-<snapatac2-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091
$ docker inspect us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091
```

## Usage

### SnapATAC2 

```bash
$ docker run --rm -it \
    -v /files:/files \
    us.gcr.io/broad-gotc-prod/samtools-bwa:1.0.0-0.7.17-1678998091 /usr/gitc/bwa mem \
    /bwa-files/ref.fa /bwa-files/reads.fq > /bwa-files/aln.sam
```

