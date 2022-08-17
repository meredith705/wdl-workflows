FROM ubuntu:20.04

MAINTAINER jmonlong@ucsc.edu

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    python3 \
    pigz \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /scripts

ADD subset_reads_by_qscore.py /scripts/

WORKDIR /home
