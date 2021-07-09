# Illumina Iaap CLI

## Quick reference

Copy and paste to pull this image

#### `us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.2-1.1.0-1625852425`

- __What is this image:__ This image is a lightweight alpine-based image for running the Illumina Iaap CLI.
- __What is the Illumina IAAP CLI:__ TODO See [here](https://emea.support.illumina.com/downloads/iaap-genotyping-cli.html) more information.
- __How to see Illumina CLI version used in image:__ Please see below.

## Versioning

The zCall image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:<image-version>-<iaap-cli-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.2-1.1.0-1625852425
$ docker inspect us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.2-1.1.0-1625852425
```

## Usage

```bash
$ docker run --rm -it \
    -v /zcall-files:/zcall-files \
    us.gcr.io/broad-gotc-prod/zcall:4.0.1-1.3-1625831535 /usr/gitc/zcall/zCall.py \
    -B /zcall-files/my.bpm.csv -G /zcall-files/my.gtc -T /zcall-files/my.thresholds.txt > /my.new.ped
```