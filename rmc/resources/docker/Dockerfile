FROM 'google/cloud-sdk:slim'
FROM 'hailgenetics/python-dill:3.11-slim'

# TODO: Double check if all of these packages are necessary
# Install packages
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
RUN wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | tee /etc/apt/keyrings/adoptium.asc && echo "deb [signed-by=/etc/apt/keyrings/adoptium.asc] https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list && apt-get update && apt install -y --no-install-recommends temurin-8-jdk

# Intall Hail
ENV HAIL_VERSION="0.2.120"
RUN python3 --version
RUN python3 -m pip install hail==${HAIL_VERSION}

# Install GCS Connector
RUN curl -sSL broad.io/install-gcs-connector > install-gcs-connector.py && \
    python3 install-gcs-connector.py --key-file-path /gsa-key/key.json

WORKDIR /home
