FROM ubuntu:20.04

# Use --build-arg #.#.# when building to use other version
ARG JAVA_MAJOR_VERSION=11
ARG ORBDETPY_VERSION=2.1.0

# Install orbdetpy dependencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    openjdk-${JAVA_MAJOR_VERSION}-jdk software-properties-common python3-pip python3-venv wget && \
    add-apt-repository ppa:deadsnakes/ppa

# Install orbdetpy into virtual environment and update orekit-data
RUN cd && \
    python3 -m venv env_orbdetpy && \
    . env_orbdetpy/bin/activate && \
    pip install orbdetpy==${ORBDETPY_VERSION} ipython && \
    python -c "from orbdetpy.astro_data import update_data; update_data();"

RUN cd && \
    wget -qO- https://github.com/ut-astria/orbdetpy/archive/refs/tags/${ORBDETPY_VERSION}.tar.gz | \
    tar -xvz -C /root/ && \
    mv orbdet* orbdetpy