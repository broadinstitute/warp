#!/bin/bash
# onCreateCommand.sh - Runs during container build to install system dependencies
set -e

echo "🔧 Setting up WARP development environment..."

# Update package lists
apt-get update
apt-get upgrade -y

# Install system dependencies
echo "📦 Installing system packages..."
apt-get install -y \
  curl \
  wget \
  git \
  build-essential \
  ca-certificates

# Install Java (required for womtool WDL validation)
echo "☕ Installing Java..."
apt-get install -y openjdk-17-jdk-headless

# Download and setup womtool
echo "🔨 Setting up womtool for WDL validation..."
WOMTOOL_VERSION="90"
WOMTOOL_DIR="/opt/womtool"
mkdir -p $WOMTOOL_DIR
cd $WOMTOOL_DIR

# Download womtool JAR if not already present
if [ ! -f "womtool-${WOMTOOL_VERSION}.jar" ]; then
  echo "Downloading womtool v${WOMTOOL_VERSION}..."
  wget -q "https://github.com/broadinstitute/cromwell/releases/download/${WOMTOOL_VERSION}/womtool-${WOMTOOL_VERSION}.jar" || \
  wget -q "https://repo1.maven.org/maven2/org/broadinstitute/cromwell/womtool/${WOMTOOL_VERSION}/womtool-${WOMTOOL_VERSION}.jar" || \
  echo "⚠️  Womtool download failed - you can run 'bash .devcontainer/download-womtool.sh' later"
fi

# Create symlink for easy access
ln -sf "womtool-${WOMTOOL_VERSION}.jar" womtool.jar

# Add womtool to PATH via wrapper script
cat > /usr/local/bin/womtool << 'EOF'
#!/bin/bash
java -jar /opt/womtool/womtool.jar "$@"
EOF
chmod +x /usr/local/bin/womtool

# Install Python 3 and pip
echo "🐍 Installing Python..."
apt-get install -y python3 python3-pip python3-venv

# Upgrade pip
python3 -m pip install --upgrade pip setuptools wheel

# Clean up apt cache
apt-get clean
rm -rf /var/lib/apt/lists/*

echo "✅ System setup complete!"
