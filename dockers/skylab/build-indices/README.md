# Build_indices

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1663605340`
`

- __What is this image:__ This image is a Debian-based custom image with STAR installed and pre-configured along with python scripts to build indices.
- __What is STAR:__ Spliced Transcripts Alignment to a Reference (STAR) is a fast RNA-seq read mapper, with support for splice-junction and fusion read detection. STAR aligns reads by finding the Maximal Mappable Prefix (MMP) hits between reads (or read pairs) and the genome, using a Suffix Array index, [more info here](https://github.com/alexdobin/STAR).
- __How to see tool version used in image:__ Please see below.

## Versioning

Build_indices uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/build-indices:<image-version>-<star-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1663605340
$ docker inspect us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1663605340
```

## Usage

### Build_indices 

```bash
$ docker run --rm -it \
    us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1663605340 \
    build-indices bash
```

Then you can exec into the container and use STAR or any of the scripts accordingly. Alternatively, you can run one-off commands by passing the command as a docker run parameter.