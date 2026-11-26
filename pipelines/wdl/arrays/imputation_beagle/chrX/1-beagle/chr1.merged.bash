#!/bin/bash
set -eu

echo    "command:  $0 $@"       1>&2
echo    "nproc:    $(nproc)"     1>&2
echo -e "date:     $(date)\n"    1>&2

echo -e "\nCopying input files:"  1>&2
time gcloud storage cp "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/hail_in_wdl_testing/outputs/aou-anvil.final.step5b.chr1.vcf.bgz" .
time gcloud storage cp "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/jar/beagle.jar" .
time gcloud storage cp "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/maps/plink.chrchr1.GRCh38.map" .

echo -e "\nRunning beagle:"       1>&2
time java -ea -Xmx780g -jar beagle.jar gt=aou-anvil.final.step5b.chr1.vcf.bgz map=plink.chrchr1.GRCh38.map out=chr1.aou-anvil window-markers=3200000

echo -e "\nCopying output files:" 1>&2
time gcloud storage cp chr1.aou-anvil.log    gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/merged.phased
time gcloud storage cp chr1.aou-anvil.vcf.gz gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/merged.phased

