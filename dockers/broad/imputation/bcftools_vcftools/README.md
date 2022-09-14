# Imputation BCFtools VCFtools

## Quick reference

Copy and paste to pull this image

#### `us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.6-1.10.2-0.1.16-1661278921`

- __What is this image:__ This image is a lightweight alpine-based image for running BCFtools and VCFtools for the [Imputation pipeline](../../../../pipelines/broad/arrays/imputation/Imputation.wdl).
- __What are BFCtools and VCFtools:__ BCFtools and VCFtools are a suite of tools for variant calling and manipulating BCFs and VCFs. See [here](https://github.com/samtools/vcftools) and [here](https://vcftools.github.io/man_latest.html) more information.
- __How to see tool versions used in image:__ Please see below.

## Versioning

The Imputation BCFtools VCFtools image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/samtools:<image-version>-<bcftools-version>-<vcftools-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.6-1.10.2-0.1.16-1661278921
$ docker inspect us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.6-1.10.2-0.1.16-1661278921
```

## Usage

### Display BCFtools default menu

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.6-1.10.2-0.1.16-1661278921 bcftools
```

### Display VCFtools default menu

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.6-1.10.2-0.1.16-1661278921 vcftools
```
