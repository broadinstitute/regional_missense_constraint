FROM 'google/cloud-sdk:slim'
FROM 'python:3.11-slim'

# install packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    curl \
    g++ \
    gawk \
    less \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    liblz4-dev \
    man-db \
    pkg-config \
    python3-venv \
    software-properties-common \
    unzip \
    wget \
    zlib1g-dev \
    && \
    # clean up apt cache
    rm -rf /var/lib/apt/lists/*

# Install Java 8 for hail
RUN mkdir -p /etc/apt/keyrings
RUN wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | tee /etc/apt/keyrings/adoptium.asc \
    && echo "deb [signed-by=/etc/apt/keyrings/adoptium.asc] https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" \
    | tee /etc/apt/sources.list.d/adoptium.list && apt-get update && apt install -y --no-install-recommends temurin-8-jdk

RUN mkdir tools
WORKDIR /tools

# Install python packages
RUN apt-get update \
    && apt-get dist-upgrade -y \
    && apt-get install -y --no-install-recommends\
    openssl \
    python3-smbus \
    python3 python3-pip

# Upgrade pip to latest version
RUN python3 -m pip install --upgrade pip

# Install hail and other python libraries
ENV HAIL_VERSION="0.2.120"
RUN python3 --version
RUN python3 -m pip install \
    wheel \
    pypandoc \
    hail==${HAIL_VERSION} \
    scipy \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    ipython \
    pybedtools \
    dill

# Install GCS Connector
RUN curl -sSL broad.io/install-gcs-connector > install-gcs-connector.py && \
    python3 install-gcs-connector.py --key-file-path /gsa-key/key.json

WORKDIR /home
