FROM centos:centos6

ENV TERM=xterm-256color

WORKDIR /usr/gitc

RUN rpm -Uvh http://download.fedoraproject.org/pub/epel/6/i386/epel-release-6-8.noarch.rpm \
  && yum update -y \
  && yum install -y curl which wget tar glib2-devel cairo-devel libtiff-devel libjpeg-devel giflib-devel libexif-devel libX11 libX11-devel libXt libXt-devel java-1.8.0-openjdk R \
  && yum groupinstall -y 'Development Tools' \
  && yum install -y zlib-devel openssl-devel sqlite-devel bzip2-devel

# put python 2.7.10 on this image so gcloud will work (deprecated as of January 2017)
RUN curl -O https://www.python.org/ftp/python/2.7.10/Python-2.7.10.tgz \
    && tar xvfz Python-2.7.10.tgz \
    && cd Python-2.7.10 \
    && ./configure --enable-shared --prefix=/usr/local LDFLAGS=-Wl,--rpath=/usr/local/lib \
    && make \
    && make altinstall \
    && cd .. \
    && rm -rf Python-2.7.10 \
    && rm Python-2.7.10.tgz

# tell cloudsdk which python to use since we dont want to mess with the system python
ENV CLOUDSDK_PYTHON=/usr/local/bin/python2.7

# download the gsutil install script
RUN curl https://dl.google.com/dl/cloudsdk/channels/rapid/install_google_cloud_sdk.bash > gsutil_install.sh

# install gsutil
RUN /bin/bash gsutil_install.sh --disable-prompts


# restart shell
RUN exec -l /bin/bash

COPY . .
