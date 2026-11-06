#!/bin/bash
set -eu

echo    "command:  $0 $@"       1>&2
echo    "nproc:    $(nproc)"     1>&2
echo -e "date:     $(date)\n"    1>&2

echo -e "\nCopying input file:"   1>&2
time gcloud storage cp "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/merged.phased/chr1.aou-anvil.vcf.gz" .

echo -e "\nRunning tabix:"        1>&2
time tabix "chr1.aou-anvil.vcf.gz"

echo -e "\nCopying output file:" 1>&2
time gcloud storage cp "chr1.aou-anvil.vcf.gz.tbi" "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/merged.phased"

