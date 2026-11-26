#!/bin/bash
set -eu

echo    "command:  $0 $@"       1>&2
echo    "nproc:    $(nproc)"     1>&2
echo -e "date:     $(date)\n"    1>&2

echo -e "\nCopying input file:" 1>&2
time gcloud storage cp "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/merged.phased/chr1.aou-anvil.vcf.gz" . 1>&2

echo -e "\nExtracting markers:"      1>&2
time zcat "chr1.aou-anvil.vcf.gz" | grep -v '^##contig' | cut -f1-9 | bgzip > chr1.aou-anvil.markers.gz

echo -e "\nCopying output file:" 1>&2
time gcloud storage cp "chr1.aou-anvil.markers.gz" "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/merged.phased" 1>&2

