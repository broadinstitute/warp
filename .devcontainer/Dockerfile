# .devcontainer/Dockerfile - Custom development container for WARP
FROM mcr.microsoft.com/devcontainers/base:ubuntu-22.04

LABEL maintainer="WARP Development"
LABEL description="Development container for WARP WDL pipelines"

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV WARP_HOME=/workspace

# Update and install system dependencies
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    curl \
    wget \
    git \
    build-essential \
    ca-certificates \
    openjdk-17-jdk-headless \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Setup womtool for WDL validation
RUN mkdir -p /opt/womtool && \
    cd /opt/womtool && \
    WOMTOOL_VERSION="90" && \
    wget -q "https://github.com/broadinstitute/cromwell/releases/download/${WOMTOOL_VERSION}/womtool-${WOMTOOL_VERSION}.jar" 2>/dev/null || \
    echo "Note: womtool will be downloaded on first use" && \
    ln -sf "womtool-${WOMTOOL_VERSION}.jar" womtool.jar 2>/dev/null || true

# Create wrapper script for womtool
RUN cat > /usr/local/bin/womtool << 'EOF' && chmod +x /usr/local/bin/womtool
#!/bin/bash
java -jar /opt/womtool/womtool.jar "$@"
EOF

# Upgrade pip and install core Python packages
RUN python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel && \
    python3 -m pip install --no-cache-dir \
    miniwdl==1.9.1 \
    pytest==7.4.3 \
    pytest-xdist==3.5.0 \
    pytest-timeout==2.2.0 \
    pytest-cov==4.1.0 \
    pyyaml==6.0.1 \
    jinja2==3.1.6 \
    jsonschema==4.20.0 \
    black==23.12.0 \
    pylint==3.0.3 \
    mypy==1.8.0 \
    ruff==0.1.11 \
    ipython \
    jupyter \
    jupyterlab

# Set working directory
WORKDIR /workspace

# Display versions on container startup
RUN cat > /usr/local/bin/show-versions << 'EOF' && chmod +x /usr/local/bin/show-versions
#!/bin/bash
echo "WARP Development Environment Versions:"
echo "======================================"
echo "Java: $(java -version 2>&1 | head -n1)"
echo "Python: $(python3 --version)"
echo "Pip: $(python3 -m pip --version | cut -d' ' -f1-2)"
echo "Git: $(git --version)"
echo ""
echo "Python Packages:"
python3 -m pip list | grep -E "miniwdl|pytest|black|mypy|ruff"
EOF

# Default command
CMD ["/bin/bash"]
