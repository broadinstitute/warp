FROM us.gcr.io/broad-dsp-gcr-public/base/jre:15-debian

ENV TERM=xterm-256color

# Change working directory to /usr/gitc/
WORKDIR /usr/gitc

# Install mono
RUN apt-get update && apt-get install -y mono-complete

# Copy everything in working dir outside container to /usr/gitc
COPY . .
