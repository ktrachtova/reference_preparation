# Use an official Ubuntu as a parent image
FROM ubuntu:22.04

# Set environment variables to non-interactive to avoid prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install R and isomiRs package
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
    unzip && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod +x blat

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    r-base-dev && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip setuptools wheel

RUN pip3 install cutadapt==5.1 biopython==1.85

RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz && \
    tar xvf cd-hit-v4.8.1-2019-0228.tar.gz --gunzip && \
    cd cd-hit-v4.8.1-2019-0228 && \
    make && \
    cd cd-hit-auxtools && \
    make

RUN git clone https://github.com/soedinglab/MMseqs2.git && \
    cd MMseqs2 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .. && \
    make -j8 && \
    make install && \
    export PATH=$(pwd)/bin/:$PATH

RUN R -e "install.packages('data.table'); install.packages('dplyr'); install.packages('ggseqlogo'); install.packages('stringr'); install.packages('ggplot2')"
    
CMD ["bash"]
