ARG UBUNTU_VER=18.04
ARG CONDA_VER=latest
ARG OS_TYPE=x86_64
ARG CROMWELL_VER=74

FROM ubuntu:${UBUNTU_VER}
COPY . /gatk-workflow/

# System packages
RUN apt-get update && apt-get install -yq curl wget jq vim nano openjdk-11-jdk
SHELL ["/bin/bash", "-c"]
# Use the above args 
ARG CONDA_VER
ARG OS_TYPE
# Install miniconda to /miniconda
RUN curl -LO "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh"
RUN bash Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh -p /miniconda -b
RUN rm Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda \
&& conda init bash \
&& conda env create -f /gatk-workflow/environment.yml

ARG CROMWELL_VER
RUN curl -LO "https://github.com/broadinstitute/cromwell/releases/download/74/cromwell-${CROMWELL_VER}.jar"
