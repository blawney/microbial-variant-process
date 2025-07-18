FROM debian:bookworm-slim

RUN apt-get update && \
    apt-get install -y default-jdk curl bzip2 build-essential zlib1g-dev python3-dev python3-pip python3-venv

RUN curl -L --output-dir /usr/local/bin -O https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar

RUN cd /tmp && curl -L -O https://github.com/lh3/bwa/archive/refs/tags/v0.7.19.tar.gz && \
  tar -xzf v0.7.19.tar.gz && \
  cd bwa-0.7.19 && \
  make && \
  mv bwa /usr/local/bin

RUN python3 -m venv /py_venv && \
  /py_venv/bin/pip install --upgrade pip && \
  /py_venv/bin/pip install CrossMap

ENV PATH="/py_venv/bin:$PATH"
