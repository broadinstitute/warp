# Illumina Iaap CLI

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.2-1.1.0-1629910298`

- __What is this image:__ This image is a lightweight alpine-based image for running the Illumina Iaap CLI.
- __What is the Illumina IAAP CLI:__ Illumina Iaap CLI is a tool for reporting genotype calls and sample-level quality metrics in various formats (Illumina GTC and PED files). See [here](https://emea.support.illumina.com/downloads/iaap-genotyping-cli.html) more information.
- __How to see Illumina CLI version used in image:__ Please see below.

## Versioning

The Illumina CLI image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:<image-version>-<iaap-cli-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.2-1.1.0-1629910298
$ docker inspect us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.2-1.1.0-1629910298
```

## Usage

```bash
$ docker run --rm -it \
    -v /illumina-files:/illumina-files \
    us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.2-1.1.0-1629910298 /usr/gitc/iaap/iaap-cli/iaap-cli \
    gencall /illumina-files/bead_pool_manifest_file /illumina-files/cluster_file .
```