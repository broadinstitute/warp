#!/bin/bash
set -eu

echo    "command:  $0 $@"       1>&2
echo    "nproc:    $(nproc)"     1>&2
echo -e "date:     $(date)\n"    1>&2

echo -e "\nCopying input files:"  1>&2
time gcloud storage cp "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/merged.phased/chr1.aou-anvil.vcf.gz" .
time gcloud storage cp "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/jar/bref3.jar" .

echo -e "\nRunning bref3:"        1>&2
time java -ea -Xmx55g -jar bref3.jar chr1.aou-anvil.vcf.gz > chr1.aou-anvil.bref3


echo -e "\nCopying output files:" 1>&2
time gcloud storage cp chr1.aou-anvil.bref3   gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/merged.phased
