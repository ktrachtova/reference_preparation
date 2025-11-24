# Use an official Ubuntu 24.04 as a parent image
FROM ubuntu:24.04

# Set environment variables to non-interactive to avoid prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# ==== Package versions (Ubuntu 24.04 / noble) ====
# cutadapt         4.4-1build2       (apt)
# python3-cutadapt 4.4-1build2       (apt)
# mmseqs2          15-6f452+ds-2     (apt)
# biopython        1.86              (pip)
# ================================================

# Base system + build deps
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        gnupg2 \
        lsb-release \
        software-properties-common \
        libpng-dev \
        libssl-dev \
        liblapack-dev \
        libopenblas-dev \
        openjdk-8-jre \
        python3 \
        python3-pip \
        python3-dev \
        build-essential \
        cython3 \
        rustc \
        wget \
        git-all \
        libtool \
        automake \
        autoconf \
        pkg-config \
        libgtextutils-dev \
        g++ \
        cmake \
        gfortran \
        libcurl4-openssl-dev \
        libxml2-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libpcre2-dev \
        unzip && \
    rm -rf /var/lib/apt/lists/*

# NOTE: This blat binary is x86_64-only and will NOT run on arm64.
RUN wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod +x blat

    
# Add CRAN repo for R (apt-key is deprecated but still works; expect a warning)
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        r-base-dev && \
    rm -rf /var/lib/apt/lists/*

# Upgrade pip toolchain
# RUN python3 -m pip install --upgrade pip setuptools wheel

# Install cutadapt from distro with explicit version (Ubuntu 24.04: 4.4-1build2)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        cutadapt=4.4-1build2 \
        python3-cutadapt=4.4-1build2 && \
    rm -rf /var/lib/apt/lists/*

# Install Biopython with explicit version from PyPI (current latest: 1.86)
RUN python3 -m pip install --break-system-packages --no-cache-dir biopython==1.86

# Install R packages
RUN R -e "install.packages('data.table'); install.packages('dplyr'); install.packages('ggseqlogo'); install.packages('stringr'); install.packages('ggplot2')"

# === MMseqs2: install official SSE2 binary (generic x86_64, NO AVX2) ===
# Clean up possible conflicting mmseqs files
RUN rm -f /usr/bin/mmseqs /bin/mmseqs || true

RUN wget https://mmseqs.com/latest/mmseqs-linux-sse2.tar.gz && \
    tar xzf mmseqs-linux-sse2.tar.gz && \
    mv mmseqs /opt/mmseqs2 && \
    ln -s /opt/mmseqs2/bin/mmseqs /usr/local/bin/mmseqs && \
    rm mmseqs-linux-sse2.tar.gz && \
    mmseqs version

CMD ["bash"]
