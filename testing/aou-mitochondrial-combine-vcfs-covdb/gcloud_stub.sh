#!/bin/bash
# Stub for 'gcloud' that redirects Google Cloud Storage operations to
# fake-gcs-server for local testing.
#
# Supported subcommands (only what build_vcf_shard_mt WDL task needs):
#
#   gcloud storage cp <local_file> gs://<bucket>/<object>
#     → uploads local_file to fake-gcs-server via curl
#
#   gcloud storage cp gs://<bucket>/<object> <local_file>
#     → downloads object from fake-gcs-server via curl
#
#   gcloud storage objects describe gs://<bucket>/<object> --format=value(md5Hash)
#     → fetches object metadata from fake-gcs-server and prints md5Hash
#
# All other invocations are forwarded to the real gcloud if it exists at
# /usr/lib/google-cloud-sdk/bin/gcloud, or silently succeed (exit 0).
#
# Required environment variable:
#   STORAGE_EMULATOR_HOST  (set by miniwdl --env flag in run_test.sh)
#     e.g. http://host.docker.internal:4443

set -euo pipefail

FAKE_GCS="${STORAGE_EMULATOR_HOST:-http://host.docker.internal:4443}"
REAL_GCLOUD="/usr/lib/google-cloud-sdk/bin/gcloud"

# Helper: parse gs://bucket/object/path into bucket + object
parse_gs_url() {
    local url="$1"
    GS_BUCKET=$(echo "$url" | sed 's|gs://||' | cut -d/ -f1)
    GS_OBJECT=$(echo "$url" | sed "s|gs://${GS_BUCKET}/||")
}

# Helper: URL-encode a string (encodes / and special chars)
url_encode() {
    python3 -c "import sys, urllib.parse; print(urllib.parse.quote(sys.argv[1], safe=''))" "$1"
}

# ── Dispatch ──────────────────────────────────────────────────────────────────

SUBCMD1="${1:-}"
SUBCMD2="${2:-}"

# storage cp <src> <dest>
# Handles both uploads (src=local, dest=gs://) and downloads (src=gs://, dest=local)
if [[ "$SUBCMD1" == "storage" && "$SUBCMD2" == "cp" ]]; then
    SRC="${3:?storage cp requires a source}"
    DEST="${4:?storage cp requires a destination}"
    
    # Upload: local file to gs://bucket/object
    if [[ "$SRC" != gs://* && "$DEST" == gs://* ]]; then
        parse_gs_url "$DEST"
        OBJECT_ENC=$(url_encode "$GS_OBJECT")
        echo "[gcloud_stub] storage cp ${SRC} → gs://${GS_BUCKET}/${GS_OBJECT}" >&2
        HTTP_STATUS=$(curl -s -o /dev/null -w "%{http_code}" \
            -X POST \
            "${FAKE_GCS}/upload/storage/v1/b/${GS_BUCKET}/o?uploadType=media&name=${OBJECT_ENC}" \
            -H "Content-Type: application/octet-stream" \
            --data-binary "@${SRC}")
        if [[ "$HTTP_STATUS" != "200" ]]; then
            echo "[gcloud_stub] ERROR: upload returned HTTP ${HTTP_STATUS}" >&2
            exit 1
        fi
        echo "[gcloud_stub] upload OK (HTTP ${HTTP_STATUS})" >&2
        exit 0
    fi
    
    # Download: gs://bucket/object to local file
    if [[ "$SRC" == gs://* && "$DEST" != gs://* ]]; then
        parse_gs_url "$SRC"
        OBJECT_ENC=$(url_encode "$GS_OBJECT")
        echo "[gcloud_stub] storage cp gs://${GS_BUCKET}/${GS_OBJECT} → ${DEST}" >&2
        mkdir -p "$(dirname "$DEST")"
        HTTP_STATUS=$(curl -s -w "%{http_code}" \
            -o "$DEST" \
            "${FAKE_GCS}/storage/v1/b/${GS_BUCKET}/o/${OBJECT_ENC}?alt=media")
        if [[ "$HTTP_STATUS" != "200" ]]; then
            echo "[gcloud_stub] ERROR: download returned HTTP ${HTTP_STATUS}" >&2
            rm -f "$DEST"
            exit 1
        fi
        echo "[gcloud_stub] download OK (HTTP ${HTTP_STATUS})" >&2
        exit 0
    fi
    
    echo "[gcloud_stub] ERROR: unsupported storage cp variant: ${SRC} → ${DEST}" >&2
    exit 1
fi

# storage objects describe gs://<bucket>/<object> [--format=value(md5Hash)]
if [[ "$SUBCMD1" == "storage" && "$SUBCMD2" == "objects" && "${3:-}" == "describe" ]]; then
    GS_URL="${4:?storage objects describe requires a gs:// URL}"
    parse_gs_url "$GS_URL"
    OBJECT_ENC=$(url_encode "$GS_OBJECT")

    # Check if caller wants just the md5Hash value
    WANT_MD5=false
    for arg in "${@:5}"; do
        if [[ "$arg" == "--format=value(md5Hash)" || "$arg" == *"md5Hash"* ]]; then
            WANT_MD5=true
            break
        fi
    done

    echo "[gcloud_stub] storage objects describe gs://${GS_BUCKET}/${GS_OBJECT}" >&2
    METADATA=$(curl -sf "${FAKE_GCS}/storage/v1/b/${GS_BUCKET}/o/${OBJECT_ENC}")
    if [[ "$WANT_MD5" == "true" ]]; then
        echo "$METADATA" | python3 -c \
            "import sys, json; print(json.load(sys.stdin).get('md5Hash', ''))"
    else
        echo "$METADATA"
    fi
    exit 0
fi

# Any other gcloud command: forward to real gcloud if available, else no-op
if [[ -x "$REAL_GCLOUD" ]]; then
    exec "$REAL_GCLOUD" "$@"
fi

echo "[gcloud_stub] ignoring unhandled command: gcloud $*" >&2
exit 0
