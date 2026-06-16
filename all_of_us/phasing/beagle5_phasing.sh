#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat 1>&2 <<'EOF'
Usage:
  run_beagle.sh \
    --vcf gs://.../input.vcf.bgz \
    --map gs://.../plink.map \
    --out-prefix chr20.aou.v9.phased \
    --out-dir gs://.../phased-test \
    [--window-markers 3500000] \
    [--java-xmx 780g]
EOF
}

window_markers="3500000"
java_xmx="780g"

vcf_gs=""
map_gs=""
out_prefix=""
out_dir=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) vcf_gs="$2"; shift 2;;
    --map) map_gs="$2"; shift 2;;
    --out-prefix) out_prefix="$2"; shift 2;;
    --out-dir) out_dir="$2"; shift 2;;
    --window-markers) window_markers="$2"; shift 2;;
    --java-xmx) java_xmx="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" 1>&2; usage; exit 2;;
  esac
done

if [[ -z "${vcf_gs}" || -z "${map_gs}" || -z "${out_prefix}" || -z "${out_dir}" ]]; then
  echo "Missing required args." 1>&2
  usage
  exit 2
fi

echo -n "command:  " 1>&2; printf '%q ' "$0" "$@" 1>&2; echo 1>&2
echo    "nproc:    $(nproc)"     1>&2
echo -e "date:     $(date)\n"    1>&2

vcf_local="$(basename "${vcf_gs}")"
map_local="$(basename "${map_gs}")"

echo -e "\nCopying input files:"  1>&2
time gcloud storage cp "${vcf_gs}" "${vcf_local}"
time gcloud storage cp "${map_gs}" "${map_local}"

echo -e "\nRunning beagle:"       1>&2
java -version
time java -ea -Xmx"${java_xmx}" -jar /opt/beagle/beagle.jar \
  gt="${vcf_local}" \
  map="${map_local}" \
  out="${out_prefix}" \
  window-markers="${window_markers}"

echo -e "\nCopying output files:" 1>&2
time gcloud storage cp "${out_prefix}.log"    "${out_dir}/"
time gcloud storage cp "${out_prefix}.vcf.gz" "${out_dir}/"
