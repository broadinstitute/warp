FROM mambaorg/micromamba:0.23.0
COPY --chown=$MAMBA_USER:$MAMBA_USER hisat3n_env.yml /tmp/hisat3n_env.yml
RUN micromamba install -y -f /tmp/hisat3n_env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN yap --version
RUN allcools --version

USER root
# default argument when not provided in the --build-arg
# to build the image with gcp, use
# docker build --build-arg gcp=true -t mapping-gcp:tag .
ARG gcp
RUN if [ "$gcp" = "true" ] ; then \
        apt-get update && \
        apt-get install -y curl gnupg && \
        echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | \
        tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
        curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
        apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
        apt-get update -y && \
        apt-get install google-cloud-sdk -y;  \
      else echo 'no gcp install';  \
    fi

RUN apt-get update; \
    apt-get install; -y make; \
    apt-get install -y git-all; \
    apt-get install -y build-essential; \
    mkdir -p ~/pkg; \
    cd ~/pkg; \
    git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n; \
    cd hisat-3n; \
    git checkout -b hisat-3n origin/hisat-3n; \
    make; \
    echo 'export PATH=$HOME/pkg/hisat-3n:$PATH' >> ~/.bashrc; \
    source ~/.bashrc 
