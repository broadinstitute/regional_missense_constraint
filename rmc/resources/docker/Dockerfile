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

#install Java 8 for hail
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
ENV HAIL_VERSION="0.2.10"
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
    dill \
    gnomad

# Install GCS Connector
RUN export SPARK_HOME=$(find_spark_home.py) && \
    curl https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.0.1.jar \
         >$SPARK_HOME/jars/gcs-connector-hadoop2-2.0.1.jar && \
    mkdir -p $SPARK_HOME/conf && \
    touch $SPARK_HOME/conf/spark-defaults.conf && \
    sed -i $SPARK_HOME/conf/spark-defaults.conf \
        -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.enable.*:spark.hadoop.google.cloud.auth.service.account.enable true:' \
        -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile.*:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile /gsa-key/key.json:'

WORKDIR /home
ENV PYSPARK_SUBMIT_ARGS="--driver-memory 8g --executor-memory 8g pyspark-shell"
