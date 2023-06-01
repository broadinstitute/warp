#!/bin/bash

TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

set -e

TMPDIR=$(mktemp -t gotc-dockerXXXXXX -d)

DOCKER_VERSION="2.5.7"
DOCKER_IMAGE_TAG="$DOCKER_VERSION-$TIMESTAMP"

PICARD_PRIVATE_VERSION="1.1448"
PICARD_PUBLIC_VERSION="2.23.8"
GATK35_VERSION="3.5-0-g36282e4"
GATK4_VERSION="4.1.8.0"
SAMTOOLS_VER="1.11"
BWA_VER="0.7.15.r1140"
K8_VER="null"
BWA_POSTALT_SCRIPT_VER="null"
TABIX_VER="0.2.5_r1005"
BGZIP_VER="1.3"
SVTOOLKIT_VER="2.00-1650"

PICARD="/seq/software/picard-public/${PICARD_PUBLIC_VERSION}/picard.jar"
GATK35="/seq/software/gotc/gatk/GenomeAnalysisTK-${GATK35_VERSION}/GenomeAnalysisTK-${GATK35_VERSION}.jar"
GATK4="/seq/software/gotc/gatk/gatk4/gatk-${GATK4_VERSION}/"
TABIX="/seq/software/picard/${PICARD_PRIVATE_VERSION}/3rd_party/tabix/tabix"
BGZIP="/seq/software/gotc/3rd_party/bgzip/bgzip"
SVTOOLKIT="/seq/software/gotc/svtoolkit/svtoolkit2.00/"

scp -T vpicard05:"$PICARD $GATK35 $TABIX $BGZIP" "${TMPDIR}"/
scp -r vpicard05:"$SVTOOLKIT" ${TMPDIR}/
scp -r vpicard05:"$GATK4" ${TMPDIR}/gatk4/

cp Dockerfile "${TMPDIR}"
cd ${TMPDIR}

mv GenomeAnalysisTK-${GATK35_VERSION}.jar GATK35.jar


# Tagged with the picard release version the jars and binaries were taken from.
docker build -t us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:${DOCKER_IMAGE_TAG} \
    --build-arg PICARD_PUBLIC_VERSION \
    --build-arg GATK35_VERSION \
    --build-arg GATK4_VERSION \
    --build-arg SAMTOOLS_VER \
    --build-arg BWA_VER \
    --build-arg TABIX_VER \
    --build-arg BGZIP_VER \
    --build-arg SVTOOLKIT_VER .
docker push us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:${DOCKER_IMAGE_TAG}

# Save current version info to local file
echo -e "$DOCKER_IMAGE_TAG\t$PICARD_PUBLIC_VERSION\tn/a\t$GATK35_VERSION\tn/a\t$GATK4_VERSION\t$SAMTOOLS_VER\t$BWA_VER\t$K8_VER\t$BWA_POSTALT_SCRIPT_VER\t$TABIX_VER\t$BGZIP_VER\t$SVTOOLKIT_VER" >> ${DIR}/build_docker_version.tsv

rm -rf ${TMPDIR}
