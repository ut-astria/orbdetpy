# A Dockerfile that will allow you install any version of 
# this repository with all its dependencies as well as install 
# the screening algorithm, caspy. 
FROM ubuntu:20.04

# Use --build-arg #.#.# when building to use other version
ARG JAVA_MAJOR_VERSION=11
ARG ORBDETPY_BRANCH=114067883eb7d5164593d68e5823a739e4fbba73
ARG CASPY_BRANCH=8d20a2d86036dde9f7c853914464388c6a325011
ARG OS_CPU_TYPE=linux-x86_64
ARG ORBDETPY_GIT_LINK=https://github.com/ut-astria/orbdetpy.git
ARG OREKIT_GIT_LINK=https://github.com/ut-astria/orbdetpy/releases/download/2.1.0/orekit-data.tar.gz
ARG CASPY_GIT_LINK=https://github.com/ut-astria/caspy.git

ENV HOME=/root

# Install orbdetpy dependencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    openjdk-${JAVA_MAJOR_VERSION}-jdk software-properties-common python3-pip python3.9-venv wget \
    git maven openssh-server curl && \
    add-apt-repository ppa:deadsnakes/ppa

# Install orbdetpy into virtual environment and update orekit-data
RUN cd && \
    python3.9 -m venv env_orbdetpy && \
    . env_orbdetpy/bin/activate && \
    pip install ipython && \
    git clone ${ORBDETPY_GIT_LINK} && \
    cd orbdetpy && \
    git checkout ${ORBDETPY_BRANCH} && \
    cd orbdetpy && \
    wget -O orekit-data ${OREKIT_GIT_LINK} && \
    tar -xzvf orekit-data && \
    mvn -e -Dos.detected.classifier=${OS_CPU_TYPE} package && \
    cd ${HOME} && \
    # Prevents bdist_wheel pip error
    pip install wheel && \
    pip install -e orbdetpy/ && \
    python -c "from orbdetpy.astro_data import update_data; update_data();"

# Clone caspy
RUN cd ${HOME} && \
    git clone ${CASPY_GIT_LINK} && \
    cd caspy && \
    git checkout ${CASPY_BRANCH} && \
    . ${HOME}/env_orbdetpy/bin/activate && \
    pip install -r requirements.txt

ENTRYPOINT . ${HOME}/env_orbdetpy/bin/activate
