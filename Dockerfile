ARG UBUNTU_VER=18.04
ARG CONDA_VER=latest
ARG OS_TYPE=x86_64
ARG CROMWELL_VER=74
ARG ME=alex

FROM ubuntu:${UBUNTU_VER}

# System packages
# RUN apt-get update && apt-get install -yq curl wget jq vim nano openjdk-11-jdk

# build image with bash as the default shell
SHELL ["/bin/bash", "-c"]
# add user in the image
RUN useradd -u 1000 ${ME} --create-home --shell /bin/bash
USER ${ME}
 
ARG CONDA_VER
ARG OS_TYPE

# Install miniconda to /miniconda
RUN curl -LO "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh" \
&& bash Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh -p /miniconda -b \
&& rm Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh

# set path for miniconda
ENV PATH=/miniconda/bin:${PATH}

# update / initialize and create desired environment 
RUN conda update -y conda \
&& conda init bash \
&& conda env create -f /gatk-workflow/environment.yml

# download cromwell jar file
ARG CROMWELL_VER
RUN curl -LO "https://github.com/broadinstitute/cromwell/releases/download/74/cromwell-${CROMWELL_VER}.jar"
