FROM 'google/cloud-sdk:slim'

# install packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
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
RUN wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | apt-key add - \
    && add-apt-repository -y https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/ \
    && apt-get update && apt-get install -y adoptopenjdk-8-hotspot

RUN mkdir tools
WORKDIR /tools

# Install python packages
RUN apt-get update \
    && apt-get dist-upgrade -y \
    && apt-get install -y --no-install-recommends\
    libc6-dev \
    libffi-dev \
    libgdbm-dev \
    liblapack-dev \
    liblapack3 \
    libncursesw5-dev \
    libopenblas-base \
    libopenblas-dev \
    libsqlite3-dev \
    libssl-dev \
    openssl \
    python-smbus \
    python3 \
    python-pip

# Upgrade pip to latest version
RUN python3 -m pip install --upgrade pip

# Install hail and other python libraries
# Install hail and other python libraries
# Including packages needed for VRS annotation functionality
# Also include version for ga4gh.vrs-python to be most recent and supporting their new variant schema
# Please update the constant VRS_VERSION in vrs_annotation_batch.py accordingly when GA4GH_VRS_VERSION is changed
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
# Install GCS Connector
RUN curl -sSL broad.io/install-gcs-connector > install-gcs-connector.py && \
    python3 install-gcs-connector.py --key-file-path /gsa-key/key.json

WORKDIR /home
