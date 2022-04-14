# GATK

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/gatk:1.3.0-4.2.6.1-1649964384`
`

- __What is this image:__ This image is a lightweight alpine-based image for running GATK (it includes both GATK4 & GATK3 -- 3.5 is included for backwards compatibility with HaplotypeCalling).
- __What is GATK:__ GATK (Genome Analysis Toolkit) is the industry standard for identifying SNPs and indels in germline DNA and RNAseq data. See [here](https://gatk.broadinstitute.org/hc/en-us) for more information.
- __How to see GATK version used in image:__ Please see below.

## Versioning

The GATK image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/gatk:<image-version>-<gatk4-version>-<unix-timestamp>` 


We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following command:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/gatk:1.3.0-4.2.6.1-1649964384
$ docker inspect us.gcr.io/broad-gotc-prod/gatk:1.3.0-4.2.6.1-1649964384
```

## Usage

### GATK4

```bash
$ docker run --rm -it \
    -v /gatk-files:/gatk-files \
    us.gcr.io/broad-gotc-prod/gatk:1.3.0-4.2.6.1-1649964384  /usr/gitc/gatk4/gatk --java-options "-Xms2000m -Xmx2500m" \
    PrintReads \
    -I /gatk-files/input_bam\
    --interval-padding 500 \
    -L /gatk-files/interval_list \
    -O /gatk-files/local.sharded.bam 
```

### GATK3 

```bash
$ docker run --rm -it \
    -v /gatk-files:/gatk-files \
    us.gcr.io/broad-gotc-prod/gatk:1.3.0-4.2.6.1-1649964384  java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m -Xmx9000m \
    -jar /usr/gitc/GATK35.jar \
    -T HaplotypeCaller \
    -R /gatk-files/ref_fasta \
    -o /gatk-files/gvcf.vcf.gz \
    -I /gatk-files/local.sharded.bam \
    -L /gatk-files/interval_list 
    -ERC GVCF
```
