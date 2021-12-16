# RNA Seq fgbio

## Quick reference

Copy and paste to pull this image

#### `us.gcr.io/broad-gotc-prod/fgbio:1.0.0-1.4.0-1638817487`

- __What is this image:__ This image is a lightweight alpine based image for running the fgbio tool set within our RNA sequencing pipeline.
- __What is fgbio:__ fgbio from Fulcrum Genomics is a set of tools to analyze genomice data with a focus on Next Generation Sequencing. See [here](https://github.com/fulcrumgenomics/fgbio) for more information.
- __How to see tool version used in image:__ Please see below.

## Versioning

fgbio uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/fgbio:<image-version>-<fgbio-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/fgbio:1.0.0-1.4.0-1638817487
$ docker inspect us.gcr.io/broad-gotc-prod/fgbio:1.0.0-1.4.0-1638817487
```

## Usage


```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/fgbio:1.0.0-1.4.0-1638817487 java -jar /usr/gitc/fgbio.jar
```