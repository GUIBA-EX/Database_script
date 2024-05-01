FROM ubuntu:22.04

RUN apt update &&\
    apt install -y default-jdk &&\
    apt install -y gcc zip tree &&\
    apt install -y git wget libncurses5-dev libbz2-dev liblzma-dev curl &&\
    apt install -y ncbi-blast+ &&\
    apt install -y python-is-python3 pip &&\
    apt install -y tabix vim less

RUN mkdir -p /home/app/ &&\
    mkdir -p /home/data &&\
    mkdir -p /home/tools/fastx_toolkit

# Install tools
WORKDIR /home/tools/fastx_toolkit
RUN wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 &&\
    tar xvf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 &&\
    rm fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
ENV PATH $PATH:/home/tools/fastx_toolkit/bin

WORKDIR /home/tools/
RUN wget https://catchenlab.life.illinois.edu/stacks/source/stacks-2.64.tar.gz &&\
    tar xvf stacks-2.64.tar.gz &&\
    rm stacks-2.64.tar.gz
WORKDIR /home/tools/stacks-2.64/
RUN ./configure &&\
    make &&\
    make install

# Copy files
WORKDIR /home/app
COPY requirements.txt ./
RUN pip install -r requirements.txt
COPY ./*.py ./
COPY classification_proposal/classification.py ./
COPY *_tmp.conf ./
