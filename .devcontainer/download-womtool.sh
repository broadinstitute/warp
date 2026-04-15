#!/bin/bash
# download-womtool.sh - Helper script to download womtool JAR file
# Run this if womtool wasn't downloaded during container setup

set -e

WOMTOOL_VERSION="${1:-90}"
WOMTOOL_DIR="/opt/womtool"

echo "📥 Downloading womtool v${WOMTOOL_VERSION}..."

mkdir -p "$WOMTOOL_DIR"
cd "$WOMTOOL_DIR"

# Try GitHub releases first
if ! wget "https://github.com/broadinstitute/cromwell/releases/download/${WOMTOOL_VERSION}/womtool-${WOMTOOL_VERSION}.jar" \
        -O "womtool-${WOMTOOL_VERSION}.jar"; then
  # Fallback to Maven repository
  echo "GitHub release not found, trying Maven repository..."
  if ! wget "https://repo1.maven.org/maven2/org/broadinstitute/cromwell/womtool/${WOMTOOL_VERSION}/womtool-${WOMTOOL_VERSION}.jar" \
          -O "womtool-${WOMTOOL_VERSION}.jar"; then
    echo "❌ Failed to download womtool"
    exit 1
  fi
fi

# Update symlink
ln -sf "womtool-${WOMTOOL_VERSION}.jar" womtool.jar

echo "✅ womtool downloaded successfully!"
echo "📌 Version: $WOMTOOL_VERSION"
echo "📂 Location: $WOMTOOL_DIR/womtool-${WOMTOOL_VERSION}.jar"
echo ""
echo "Test it with: womtool validate --help"
