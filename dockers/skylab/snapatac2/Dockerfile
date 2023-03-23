# Adding a platform tag to ensure that images built on ARM-based machines
# (ex. M-series macs) won't cause issues with our automated PR test suite.
# However, this is not relevant for automated builds in a CI/CD pipeline that is AMD-based.
FROM --platform="linux/amd64" python:3.9

ENV TERM=xterm-256color \
        TINI_VERSION=v0.19.0\
        PATH=$PATH:/usr/gitc

ARG SNAPATAC2_VERSION=2.2.0

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        SNAPATAC2_VERSION=${SNAPATAC2_VERSION}

WORKDIR /usr/gitc

RUN pip install cmake; \
    curl https://sh.rustup.rs -sSf | sh -s -- -y; \
    source "$HOME/.cargo/env"; \
    pip install snapatac2; \
    # Install macs2
    pip install macs2; \
    # Install TINI
    mkdir temp; \
    cd temp; \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;


# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]