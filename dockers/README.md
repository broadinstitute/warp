# Docker Style Guide 

This style guide provides formatting guidelines and best practices for writing Dockerfiles in WARP.

## :book: Table of Contents

* [Overview](#overview)  
* [Goals](#goals)
  * [Small images](#small) 
    * [Alpine base](#alpine)
    * [Minimal RUN steps](#minimal-run)
  * [Publicly accesible](#publicly)
  * [Image scanning](#scanning)
  * [Semantic tagging](#semantic)
  * [Proper signal handling](#signal)
* [Build Scripts](#build)
* [Formatting](#formatting)

## <a name="overview"></a> Overview

WARP maintains a collection of docker images which are used as execution environements for various cloud-optimized data processing pipelines. Many of these image require specific sets of tools and dependencies to run and can be thought of as _custom_ images rather than traditional application images. 
Building and maintaining these images can be challenging; this document provides a set of guidelines to assist in developing docker iamges for WARP.

## <a name="goals"></a> Goals

The following are some goals/guidelines we want to strive for when writing our Dockerfiles.

### <a name="small"></a> Small images

Building a smaller image offers advantages such as faster upload and download times along with reduced storage costs and minimized attack vector. Two of the easiest ways to minimize the size of your image is to use a small base image and to reduce the number of layers in your image.

#### <a name="alpine"></a> Alpine base

The easiest way to have a small image is to use an [Alpine](https://alpinelinux.org/) base image. Alpine linux, compared to Debian, RHEL etc., is designed specifically for security and resource efficiency, and lends itself perfectly to be used as a building block for Docker images.

Along with being a small base, Alpine also has built in deletion of package index and provides [tini](https://github.com/krallin/tini) natively through APK.

There are some instances where a Debian base image is unavoidable, specifically in the case where dependencies don't exists in APK. It is suggested that you only go to a Debian base as a last resort.


##### :eyes: Example
```dockerfile

# GOOD
FROM alpine:3.9

RUN set eux; \
        apk add --no-cache \
            curl \
            bash \

# OKAY, NOT GREAT
FROM python:debian

RUN set eux; \
        apt-get-update; \
        apt-get install -y \
            curl \
            bash \
    ; \
# Must clean up cache manually with Debian
    apt-get clean && rm -rf /var/lib/apt/list/*
```