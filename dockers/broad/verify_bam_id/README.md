# VerifyBamID

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/verify-bam-id:1.0.1-c1cba76e979904eb69c31520a0d7f5be63c72253-1639071840`

- __What is this image:__ This image is a debian-based image for running the VerifyBamID tool.
- __What is VerifyBamID:__ VerifyBamId is a tool for DNA contamination estimation from sequence reads using ancestry-agnostic method. See [here](https://github.com/Griffan/VerifyBamID) more information.
- __How to see VerifyBamID version used in image:__ Please see below.

## Versioning

The VerifyBamID image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/verify-bam-id:<image-version>-<VerifyBamID-git-hash>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/verify-bam-id:1.0.1-c1cba76e979904eb69c31520a0d7f5be63c72253-1639071840
$ docker inspect us.gcr.io/broad-gotc-prod/verify-bam-id:1.0.1-c1cba76e979904eb69c31520a0d7f5be63c72253-1639071840
```

## Usage

```bash
$ docker run --rm -it \
    -v /bam-files:/bam-files \
    us.gcr.io/broad-gotc-prod/verify-bam-id:1.0.1-c1cba76e979904eb69c31520a0d7f5be63c72253-1639071840 /usr/gitc/VerifyBamID \
    --Verbose \
    --Output output_prefix \
    --BamFile /bam-files/input_bam \
    --Reference /bam_files/ref_fasta \
    --UDPath /bam_files/contamination_sites_ud \
    --MeanPath /bam_files/contamination_site_mu \
    --BedPath /bam_files/contamination_site_bed \
    1>/dev/null
```


