#!/bin/bash

###############################################################################
# Setup script for fake-gcs-server with test data
#
# This script:
# 1. Creates the test bucket in fake-gcs-server
# 2. Uploads coverage TSV files to the bucket
#
# Prerequisites:
# - fake-gcs-server must be running on localhost:4443
#   docker run -d -p 4443:4443 --name fake-gcs \
#     fsouza/fake-gcs-server -scheme http \
#     -external-url http://host.docker.internal:4443 \
#     -public-host host.docker.internal:4443
#
# IMPORTANT: Use host.docker.internal (not localhost) for -external-url and
# -public-host. Docker Swarm task containers cannot resolve 'localhost' as the
# Mac host — they'd connect to themselves. host.docker.internal is Docker
# Desktop's DNS alias for the host machine and works from both host and containers.
###############################################################################

set -e

GCS_HOST="http://localhost:4443"
BUCKET_NAME="fake-bucket"

# Default to current directory if no argument provided
# This assumes the script is run from within the repo or repo root
REPO_ROOT="${1:-.}"

# If REPO_ROOT is just ".", resolve it to the actual directory
# and ensure we can find testing/ subdirectory
if [ "$REPO_ROOT" = "." ]; then
    REPO_ROOT="$(cd "$REPO_ROOT" && pwd)"
fi

echo "Setting up fake-gcs-server for testing..."
echo "Repository root: $REPO_ROOT"
echo ""

# Check if fake-gcs-server is running
echo "Checking for fake-gcs-server at $GCS_HOST..."
if ! curl -s -f -o /dev/null "$GCS_HOST/storage/v1/b" 2>/dev/null; then
    echo "✗ fake-gcs-server is not running at $GCS_HOST"
    echo ""
    echo "Start it with:"
    echo "  docker run -d -p 4443:4443 --name fake-gcs \\"
    echo "    fsouza/fake-gcs-server -scheme http \\"
    echo "    -external-url http://host.docker.internal:4443 \\"
    echo "    -public-host host.docker.internal:4443"
    exit 1
fi
echo "✓ fake-gcs-server is running"
echo ""

# Create bucket
echo "Creating bucket '$BUCKET_NAME'..."
if curl -s -X POST "$GCS_HOST/storage/v1/b" \
    -H "Content-Type: application/json" \
    -d "{\"name\":\"$BUCKET_NAME\"}" 2>/dev/null | grep -q "name"; then
    echo "✓ Bucket created (or already exists)"
else
    echo "⚠ Could not create bucket (may already exist)"
fi
echo ""

# Upload coverage files
echo "Uploading coverage files..."
for coverage_file in "$REPO_ROOT/testing/coverage_s001.tsv" "$REPO_ROOT/testing/coverage_s002.tsv"; do
    if [ ! -f "$coverage_file" ]; then
        echo "✗ Missing file: $coverage_file"
        exit 1
    fi
    
    filename=$(basename "$coverage_file")
    echo "  Uploading $filename..."
    
    if curl -s -X POST "$GCS_HOST/upload/storage/v1/b/$BUCKET_NAME/o?uploadType=media&name=$filename" \
        --data-binary "@$coverage_file" > /dev/null 2>&1; then
        echo "  ✓ Uploaded $filename"
    else
        echo "  ✗ Failed to upload $filename"
        exit 1
    fi
done
echo ""

# Upload VCF files
echo "Uploading VCF files..."
for vcf_file in "$REPO_ROOT/testing/s001.vcf" "$REPO_ROOT/testing/s002.vcf"; do
    if [ ! -f "$vcf_file" ]; then
        echo "✗ Missing file: $vcf_file"
        exit 1
    fi

    filename=$(basename "$vcf_file")
    echo "  Uploading $filename..."

    if curl -s -X POST "$GCS_HOST/upload/storage/v1/b/$BUCKET_NAME/o?uploadType=media&name=$filename" \
        --data-binary "@$vcf_file" > /dev/null 2>&1; then
        echo "  ✓ Uploaded $filename"
    else
        echo "  ✗ Failed to upload $filename"
        exit 1
    fi
done
echo ""

# Create and upload MT shard tar file for merge tests
echo "Creating MT shard tar files..."
mkdir -p "$REPO_ROOT/tmp/mt_shards"

# Create a minimal MT directory structure for testing
MT_SHARD_DIR="$REPO_ROOT/tmp/mt_shards/vcf_shard_00000"
mkdir -p "$MT_SHARD_DIR"

# Create a minimal VCF file for testing
cat > "$MT_SHARD_DIR/vcf_shard_00000.mt" <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
MT	100	.	A	G	.	PASS	.
EOF

# Create tar.gz
tar -czf "$REPO_ROOT/tmp/mt_shards/vcf_shard_00000_shard_mt.tar.gz" -C "$REPO_ROOT/tmp/mt_shards" vcf_shard_00000/

if [ -f "$REPO_ROOT/tmp/mt_shards/vcf_shard_00000_shard_mt.tar.gz" ]; then
    echo "  ✓ Created vcf_shard_00000_shard_mt.tar.gz"
    
    # Upload to fake bucket
    echo "  Uploading vcf_shard_00000_shard_mt.tar.gz to bucket..."
    if curl -s -X POST "$GCS_HOST/upload/storage/v1/b/$BUCKET_NAME/o?uploadType=media&name=vcf_shard_00000_shard_mt.tar.gz" \
        --data-binary "@$REPO_ROOT/tmp/mt_shards/vcf_shard_00000_shard_mt.tar.gz" > /dev/null 2>&1; then
        echo "  ✓ Uploaded vcf_shard_00000_shard_mt.tar.gz"
    else
        echo "  ✗ Failed to upload vcf_shard_00000_shard_mt.tar.gz"
        exit 1
    fi
else
    echo "  ✗ Failed to create vcf_shard_00000_shard_mt.tar.gz"
    exit 1
fi
echo ""

# Upload artifact-prone-sites BED file
echo "Uploading artifact-prone-sites BED file..."
blacklist_file="$REPO_ROOT/testing/blacklist_sites.hg38.chrM.bed"
if [ ! -f "$blacklist_file" ]; then
    echo "✗ Missing file: $blacklist_file"
    exit 1
fi

filename="blacklist_sites.hg38.chrM.bed"
echo "  Uploading $filename..."

if curl -s -X POST "$GCS_HOST/upload/storage/v1/b/$BUCKET_NAME/o?uploadType=media&name=$filename" \
    --data-binary "@$blacklist_file" > /dev/null 2>&1; then
    echo "  ✓ Uploaded $filename"
else
    echo "  ✗ Failed to upload $filename"
    exit 1
fi
echo ""

# Upload mock resource files for add_annotations.py
echo "Uploading mock resource files for annotations..."

# Function to upload a file with a specific GCS object path
upload_file_to_gcs() {
    local local_path="$1"
    local gcs_path="$2"  # e.g., "resources/mitochondria/phylotree/rCRS-centered_phylo_vars_final_update.txt"
    
    if [ ! -f "$local_path" ]; then
        echo "✗ Missing file: $local_path"
        return 1
    fi
    
    echo "  Uploading to gs://$BUCKET_NAME/$gcs_path..."
    
    # URL encode the GCS path properly (including slashes as %2F)
    # Use Python for reliable URL encoding since bash has limitations
    gcs_path_encoded=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$gcs_path', safe=''))")
    
    # Build the upload URL and show it for debugging
    upload_url="$GCS_HOST/upload/storage/v1/b/$BUCKET_NAME/o?uploadType=media&name=$gcs_path_encoded"
    # echo "    DEBUG: URL = $upload_url" >&2
    # echo "    DEBUG: Encoded path = $gcs_path_encoded" >&2
    
    if curl -s -X POST "$upload_url" \
        --data-binary "@$local_path" > /dev/null 2>&1; then
        echo "  ✓ Uploaded $gcs_path"
        return 0
    else
        echo "  ✗ Failed to upload $gcs_path"
        # Try again with verbose output to see the actual error
        echo "    Retrying with verbose output..."
        curl -v -X POST "$upload_url" --data-binary "@$local_path" 2>&1 | head -20
        return 1
    fi
}

# Upload phylotree mock
upload_file_to_gcs "$REPO_ROOT/testing/phylotree_mock.tsv" "resources/mitochondria/phylotree/rCRS-centered_phylo_vars_final_update.txt" || exit 1

# Upload variant context mock
upload_file_to_gcs "$REPO_ROOT/testing/variant_context_mock.tsv" "resources/mitochondria/variant_context/chrM_pos_ref_alt_context_categories.txt" || exit 1

# Upload PON-mt-tRNA mock
upload_file_to_gcs "$REPO_ROOT/testing/pon_mt_trna_mock.tsv" "resources/mitochondria/trna_predictions/pon_mt_trna_predictions_08_27_2020.txt" || exit 1

# Upload MitoTIP mock
upload_file_to_gcs "$REPO_ROOT/testing/mitotip_mock.tsv" "resources/mitochondria/trna_predictions/mitotip_scores_08_27_2020.txt" || exit 1

echo ""

echo "✓ Setup complete!"
echo ""
echo "You can now run:"
echo "  ~/codes/misc/warp/run_test.sh ~/codes/misc/warp"
