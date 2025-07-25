FROM debian:bookworm-slim

RUN apt-get update && \
    apt-get install -y \
    default-jdk \
    curl \
    bzip2 \
    build-essential \
    zlib1g-dev \
    python3-dev \
    python3-pip \
    python3-venv

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

ENV MINIFORGE_VERSION=23.3.1-0
ENV CONDA_DIR=/opt/conda
RUN curl -fsSL https://github.com/conda-forge/miniforge/releases/download/$MINIFORGE_VERSION/Miniforge3-$MINIFORGE_VERSION-$(uname)-$(uname -m).sh -o /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p $CONDA_DIR && \
    rm /tmp/miniforge.sh

# Activate conda by default and add to PATH
ENV PATH=$CONDA_DIR/bin:$PATH

# Ensure conda and mamba are up-to-date
RUN conda update -y conda && \
    conda install -y mamba -c conda-forge

# Install packages with mamba
RUN mamba install -y -c bioconda -c conda-forge pbmm2