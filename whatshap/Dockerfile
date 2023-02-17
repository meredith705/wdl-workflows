FROM ubuntu:20.04

MAINTAINER jmonlong@ucsc.edu

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential \
    python3-dev \
    python3-pip \
    bzip2 \
    libbz2-dev \
    liblzma-dev \
    bcftools \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install whatshap==1.7

ADD chr_lengths /opt/whatshap/

ENV PATH=$PATH:/root/.local/bin
