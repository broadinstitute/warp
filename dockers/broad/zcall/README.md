# zCall

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/zcall:4.0.1-1.3-1629910423`

- __What is this image:__ This image is a lightweight alpine-based image for running Zcall.
- __What is zCall:__ zCall is a variant caller specifically designed for calling rare single nucleotide polymorphisms (SNPs) from array-based technology. See [here](https://github.com/jigold/zCall) more information.
- __How to see zCall version used in image:__ Please see below.

## Versioning

The zCall image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/zcall:<image-version>-<zcall-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/zcall:4.0.1-1.3-1629910423
$ docker inspect us.gcr.io/broad-gotc-prod/zcall:4.0.1-1.3-1629910423
```

## Usage

```bash
$ docker run --rm -it \
    -v /zcall-files:/zcall-files \
    us.gcr.io/broad-gotc-prod/zcall:4.0.1-1.3-1629910423 python /usr/gitc/zcall/zCall.py \
    -B /zcall-files/my.bpm.csv -G /zcall-files/my.gtc -T /zcall-files/my.thresholds.txt > /my.new.ped
```