#!/bin/bash

###############################################################################
# IMPORTANT: This script is designed to run on the HOST MACHINE (macOS)
# and NOT from within the dev container.
#
# The script:
# - Runs miniwdl on the host to avoid docker-in-docker complexity
# - Manages Docker configuration for proper container access
# - Uses the host's Docker daemon to execute workflow tasks
#
# DO NOT run this from the dev container. Instead, run it from your local
# machine where miniwdl and Docker are installed.
###############################################################################

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
REPO_DIR="${1:-.}"
WORKFLOW="${2:-all_of_us/mitochondria/mt_coverage_merge.wdl}"
INPUTS="${3:-testing/inputs.json}"
SKIP_CLEANUP="${SKIP_CLEANUP:-false}"
# RERUN is a legacy alias for SKIP_CLEANUP; the call cache (always enabled)
# handles skipping completed tasks automatically across every run.
if [ "${RERUN:-false}" == "true" ]; then
    SKIP_CLEANUP=true
fi

# Control flags for Docker and fake-gcs setup
SKIP_BUILD="${SKIP_BUILD:-false}"
RESET_FAKE_GCS="${RESET_FAKE_GCS:-false}"

# Docker image configuration
# The local test image extends the production Hail image with:
#   - A VCF pre-localizer wrapper (downloads gs://fake-bucket/* VCFs from fake-gcs)
#   - A gcloud stub (redirects storage cp / objects describe to fake-gcs via curl)
# This lets build_vcf_shard_mt run locally without real GCP credentials.
HAIL_PROD_IMAGE="us.gcr.io/broad-gotc-prod/aou-mitochondrial-combine-vcfs-covdb:1.0.1"
HAIL_TEST_IMAGE="aou-mitochondrial-combine-vcfs-covdb:local-test"
ORIG_HAIL_IMAGE_ID=""  # populated by shadow_prod_image; used by restore_prod_image

# Functions
cleanup_old_runs() {
    echo -e "${YELLOW}Cleaning up old miniwdl runs...${NC}"
    # miniwdl names run dirs YYYYMMDD_HHMMSS_<workflow>
    rm -rf [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_*_mt_coverage_merge 2>/dev/null || true
    echo -e "${GREEN}Cleanup done${NC}"
}

backup_docker_config() {
    if [ -f ~/.docker/config.json ]; then
        echo -e "${YELLOW}Backing up Docker config...${NC}"
        mv ~/.docker/config.json ~/.docker/config.json.bak
    fi
}

restore_docker_config() {
    if [ -f ~/.docker/config.json.bak ]; then
        echo -e "${YELLOW}Restoring Docker config...${NC}"
        mv ~/.docker/config.json.bak ~/.docker/config.json
    fi
}

show_error() {
    local stderr_file="$1"
    if [ -f "$stderr_file" ]; then
        echo -e "${RED}━━━━━━━━━━━━━━━━━━━━ ERROR OUTPUT ━━━━━━━━━━━━━━━━━━━━${NC}"
        cat "$stderr_file"
        echo -e "${RED}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    fi
}

find_latest_run_dir() {
    ls -dt [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_*_mt_coverage_merge 2>/dev/null | head -1
}

# Patch WDL file to use test docker images instead of production images
# This creates a temporary copy of the WDL with docker tags replaced.
# Only patches docker images that have test versions available (indicated by
# the presence of a corresponding test image directory in testing/).
patch_wdl_for_testing() {
    local original_wdl="$1"
    local patched_wdl="$REPO_DIR/tmp/$(basename "$original_wdl" .wdl).test.wdl"
    
    if [ ! -f "$original_wdl" ]; then
        echo -e "${RED}✗ WDL file not found: $original_wdl${NC}"
        return 1
    fi
    
    mkdir -p "$REPO_DIR/tmp"
    
    # Copy and patch docker image references
    cp "$original_wdl" "$patched_wdl"
    
    # Replace production docker images with test images only if a test version exists.
    # Check for corresponding test image directories in testing/ to determine which
    # images have test versions available.
    if [ -d "$REPO_DIR/testing/aou-mitochondrial-combine-vcfs-covdb" ]; then
        sed -i '' 's|us\.gcr\.io/broad-gotc-prod/aou-mitochondrial-combine-vcfs-covdb:[0-9.]*|aou-mitochondrial-combine-vcfs-covdb:local-test|g' "$patched_wdl"
    fi
    
    echo "$patched_wdl"
}

# Ensure required test objects exist in fake-gcs-server.
# If any are missing, run setup_fake_gcs.sh to upload them.
# If RESET_FAKE_GCS is set, always re-run setup to overwrite existing data.
ensure_fake_gcs_data() {
    if [ "$RESET_FAKE_GCS" == "true" ]; then
        echo -e "${YELLOW}Forcing reset of fake-gcs data (RESET_FAKE_GCS=true)...${NC}"
        "$REPO_DIR/testing/setup_fake_gcs.sh" "$REPO_DIR"
        return 0
    fi
    
    local objects=("coverage_s001.tsv" "coverage_s002.tsv" "s001.vcf" "s002.vcf" "vcf_shard_00000_shard_mt.tar.gz" "blacklist_sites.hg38.chrM.bed")
    local needs_setup=false
    for obj in "${objects[@]}"; do
        if ! curl -sf "http://localhost:4443/storage/v1/b/fake-bucket/o/${obj}" > /dev/null 2>&1; then
            needs_setup=true
            break
        fi
    done
    if [ "$needs_setup" == "true" ]; then
        echo -e "${YELLOW}Test data missing from fake-gcs — running setup...${NC}"
        "$REPO_DIR/testing/setup_fake_gcs.sh" "$REPO_DIR"
    else
        echo -e "${GREEN}✓ Test data already in fake-gcs${NC}"
    fi
}

# Build the local test Docker image for build_vcf_shard_mt.
# The image extends the production Hail image and adds a VCF pre-localizer
# wrapper and a gcloud stub so the task works against fake-gcs-server.
#
# Both files are inlined in the Dockerfile so Docker's normal layer cache is
# correct: a layer is reused only when its RUN command string is unchanged,
# i.e. only when the wrapper/stub content is unchanged.  This keeps the image
# digest stable run-to-run, which is required for miniwdl's call cache to hit.
build_test_image() {
    echo -e "${YELLOW}Building local test image (${HAIL_TEST_IMAGE})...${NC}"
    if docker build -t "$HAIL_TEST_IMAGE" "$REPO_DIR/testing/aou-mitochondrial-combine-vcfs-covdb" 2>&1 \
            | tail -5; then
        echo -e "${GREEN}✓ Test image built${NC}"
    else
        echo -e "${RED}✗ Failed to build test image${NC}"
        return 1
    fi
}

# Replace the production Hail image with the local test image by removing the
# old tag entirely and retagging with the production name. This ensures miniwdl
# uses the test image regardless of whether it uses tags or digest-based lookup.
shadow_prod_image() {
    # First, save the original image ID if the tag exists
    ORIG_HAIL_IMAGE_ID=$(docker inspect --format='{{.Id}}' "$HAIL_PROD_IMAGE" 2>/dev/null || echo "")
    
    # Completely remove the old tag so there's no ambiguity
    if [ -n "$ORIG_HAIL_IMAGE_ID" ]; then
        echo -e "${YELLOW}Removing old ${HAIL_PROD_IMAGE}...${NC}"
        docker rmi "$HAIL_PROD_IMAGE" 2>/dev/null || true
    fi
    
    # Now tag our test image as the production image (only tag, no alias)
    echo -e "${YELLOW}Tagging test image as ${HAIL_PROD_IMAGE}...${NC}"
    docker tag "$HAIL_TEST_IMAGE" "$HAIL_PROD_IMAGE" 2>/dev/null
    echo -e "${GREEN}→ Production image replaced with test image${NC}"
}

# Restore the original production Hail image by:
# 1. Removing the test image tag
# 2. Re-pulling or re-tagging the original if we saved it
restore_prod_image() {
    echo -e "${YELLOW}Restoring Docker images...${NC}"
    
    # Remove the test image tag (production name) so it doesn't shadow
    docker rmi "$HAIL_PROD_IMAGE" 2>/dev/null || true
    
    # If we had an original, we'd need to re-pull it (docker doesn't store old image
    # IDs well after removal, so for simplicity we just remove and let next run fetch)
    if [ -n "$ORIG_HAIL_IMAGE_ID" ]; then
        echo -e "${YELLOW}Original image was replaced; you may need to re-pull on next run${NC}"
    fi
    echo -e "${GREEN}→ Production image tag cleaned up${NC}"
}

# Main
main() {
    cd "$REPO_DIR"
    
    # Check for fake-gcs-server
    echo -e "${YELLOW}Checking for fake-gcs-server at localhost:4443...${NC}"
    if ! curl -s -f -o /dev/null http://localhost:4443/storage/v1/b 2>/dev/null; then
        echo -e "${RED}✗ fake-gcs-server is not running!${NC}"
        echo -e "${YELLOW}Start it with:${NC}"
        echo "  docker run -d -p 4443:4443 --name fake-gcs fsouza/fake-gcs-server -scheme http -external-url http://host.docker.internal:4443 -public-host host.docker.internal:4443"
        echo ""
        restore_docker_config
        return 1
    fi
    echo -e "${GREEN}✓ fake-gcs-server is running${NC}"
    echo ""

    # Upload test data if missing from fake-gcs
    ensure_fake_gcs_data
    echo ""

    # Optional cleanup (skipped when SKIP_CLEANUP=true or RERUN=true)
    if [ "$SKIP_CLEANUP" != "true" ]; then
        cleanup_old_runs
    fi
    
    # Setup
    backup_docker_config
    
    # Build Docker image (unless SKIP_BUILD=true)
    if [ "$SKIP_BUILD" != "true" ]; then
        build_test_image || { restore_docker_config; return 1; }
    else
        echo -e "${YELLOW}Skipping Docker image build (SKIP_BUILD=true)${NC}"
        echo -e "${YELLOW}Using existing image ${HAIL_TEST_IMAGE}${NC}"
    fi
    
    shadow_prod_image
    
    # Run
    echo -e "${YELLOW}Running miniwdl workflow...${NC}"
    echo "  Workflow: $WORKFLOW"
    echo "  Inputs: $INPUTS"
    echo ""
    
    # Configure Docker and storage for host machine execution (macOS)
    # STORAGE_EMULATOR_HOST points to fake-gcs-server which mocks Google Cloud Storage
    # On Docker Desktop for Mac, host.docker.internal resolves to the host machine
    # This allows containers to reach services running on localhost:4443 of the host
    export STORAGE_EMULATOR_HOST=http://host.docker.internal:4443
    export DOCKER_HOST=unix://$HOME/.docker/run/docker.sock
    
    # --env passes the variable INTO each Docker container's environment.
    # Without this, STORAGE_EMULATOR_HOST is only in the host shell and
    # gcsfs inside mt-coverage-db:1.0.0 will connect to real GCS and get a 401.

    # Write a miniwdl cfg enabling the call cache.
    # The cache is keyed by input hash + image digest; tasks whose inputs and
    # image are unchanged are retrieved from cache and not re-executed.
    # This applies on every run — no special flag needed to resume.
    MINIWDL_CFG="$REPO_DIR/tmp/miniwdl.cfg"
    CALL_CACHE_DIR="$REPO_DIR/tmp/miniwdl-call-cache"
    mkdir -p "$CALL_CACHE_DIR"
    cat > "$MINIWDL_CFG" << CFGEOF
[call_cache]
get = true
put = true
dir = $CALL_CACHE_DIR
CFGEOF

    # Patch WDL to use test docker images
    echo -e "${YELLOW}Patching WDL docker images for testing...${NC}"
    PATCHED_WORKFLOW=$(patch_wdl_for_testing "$WORKFLOW") || { restore_docker_config; return 1; }
    echo -e "${GREEN}✓ WDL patched: $PATCHED_WORKFLOW${NC}"
    echo ""

    if miniwdl run "$PATCHED_WORKFLOW" -i "$INPUTS" \
        --cfg "$MINIWDL_CFG" \
        --env STORAGE_EMULATOR_HOST=http://host.docker.internal:4443; then
        echo -e "${GREEN}✓ Workflow completed successfully!${NC}"
        restore_prod_image
        restore_docker_config
        return 0
    else
        EXIT_CODE=$?
        echo -e "${RED}✗ Workflow failed with exit code $EXIT_CODE${NC}"
        
        # Find the latest run and show errors from the FAILED task.
        # Use the most recently modified stderr.txt (the failed task's, not a
        # succeeded task that happens to sort first alphabetically).
        LATEST_RUN=$(find_latest_run_dir)
        if [ -n "$LATEST_RUN" ]; then
            FAILED_STDERR=$(find "$LATEST_RUN" -name "stderr.txt" -printf "%T@ %p\n" 2>/dev/null \
                | sort -rn | head -1 | cut -d' ' -f2-)
            # macOS find lacks -printf; fall back to ls -t
            if [ -z "$FAILED_STDERR" ]; then
                FAILED_STDERR=$(find "$LATEST_RUN" -name "stderr.txt" 2>/dev/null \
                    | xargs ls -t 2>/dev/null | head -1)
            fi
            if [ -n "$FAILED_STDERR" ]; then
                echo ""
                show_error "$FAILED_STDERR"
            fi
        fi

        restore_prod_image
        restore_docker_config
        return $EXIT_CODE
    fi
}

# Help
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    cat << EOF
Usage: $0 [REPO_DIR] [WORKFLOW] [INPUTS]

Run miniwdl workflow on the HOST MACHINE (not from dev container).

PREREQUISITES
=============
Before running this script, ensure fake-gcs-server is running:

  docker run -d -p 4443:4443 --name fake-gcs \\
    fsouza/fake-gcs-server -scheme http \\
    -external-url http://host.docker.internal:4443 \\
    -public-host host.docker.internal:4443

This script:
  - Runs miniwdl from your local machine to avoid docker-in-docker complexity
  - Automatically manages Docker config for proper container access
  - Uses the host's Docker daemon to execute workflow tasks
  - Connects to fake-gcs-server to mock Google Cloud Storage for testing
  - Cleans up old run directories and restores Docker config on completion

Arguments:
  REPO_DIR   Repository root directory (default: current directory)
  WORKFLOW   Path to WDL workflow (default: all_of_us/mitochondria/mt_coverage_merge.wdl)
  INPUTS     Path to inputs JSON (default: testing/inputs.json)

Environment Variables:
  SKIP_CLEANUP         Set to "true" to preserve old run directories (useful for
                       inspecting outputs). By default, old run dirs are deleted
                       before each run.
  
  SKIP_BUILD           Set to "true" to skip rebuilding the Docker test image.
                       Use this if the image hasn't changed and you want a faster
                       run. Requires the image to already exist.
  
  RESET_FAKE_GCS       Set to "true" to force re-run of setup_fake_gcs.sh after
                       Docker rebuild, overwriting any existing test data in the
                       fake-gcs-server. Useful when wrappers change or test data
                       fixtures are updated. By default, setup only runs if data
                       is missing.
  
  RERUN                Alias for SKIP_CLEANUP=true (legacy name, kept for convenience).

  Note: miniwdl's call cache (in tmp/miniwdl-call-cache/) is always active.
  Tasks whose inputs and image are unchanged are skipped automatically on every
  run — no special flag is needed to resume after a failure.

Configuration
==============
The script uses STORAGE_EMULATOR_HOST=http://host.docker.internal:4443 to direct GCS 
requests from Docker containers to fake-gcs-server on the host machine. 
(host.docker.internal is Docker Desktop for Mac's special DNS name for the host machine)

This allows testing workflows that reference GCS paths (e.g., gs://fake-bucket/coverage_s001.tsv)
without requiring actual Google Cloud credentials. Test data is uploaded automatically
if not already present in fake-gcs-server (no need to run setup_fake_gcs.sh manually).

Examples:
  # Normal run (fresh image, auto-uploads data, cleans old run dirs)
  $0 ~/codes/misc/warp

  # Retry after a failure — completed tasks are skipped via call cache
  $0 ~/codes/misc/warp

  # Keep old run dirs for inspection (call cache still skips completed tasks)
  SKIP_CLEANUP=true $0 ~/codes/misc/warp

  # Fast run without rebuilding the Docker image (if unchanged)
  SKIP_BUILD=true $0 ~/codes/misc/warp

  # Force re-setup of fake-gcs data (overwrites existing, useful for wrapper changes)
  RESET_FAKE_GCS=true $0 ~/codes/misc/warp

  # Combine options: skip build + reset fake-gcs + keep old runs
  SKIP_BUILD=true RESET_FAKE_GCS=true SKIP_CLEANUP=true $0 ~/codes/misc/warp

Notes:
  - This must be executed from your local machine, not from the dev container
  - fake-gcs-server must be running before starting the workflow
  - The script will check for fake-gcs-server availability at startup
  - Test data (VCFs, coverage TSVs) is auto-uploaded to fake-gcs if missing

EOF
    exit 0
fi

main
