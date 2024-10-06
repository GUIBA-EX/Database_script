FROM ubuntu:22.04

# 添加架构检测
ARG TARGETARCH

RUN apt update && apt install -y \
    default-jdk \
    gcc \
    zip \
    tree \
    git \
    wget \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    curl \
    ncbi-blast+ \
    python-is-python3 \
    pip \
    tabix \
    vim \
    less

RUN mkdir -p /home/app/ /home/data /home/tools

# 安装 stacks
WORKDIR /home/tools/
RUN wget https://catchenlab.life.illinois.edu/stacks/source/stacks-2.64.tar.gz && \
    tar xvf stacks-2.64.tar.gz && \
    rm stacks-2.64.tar.gz && \
    cd stacks-2.64 && \
    ./configure && \
    make && \
    make install

# 安装 fastp
WORKDIR /home/tools/
RUN wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp && \
    mv ./fastp /usr/local/bin/

# 安装 seqkit
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "aarch64" ]; then \
        SEQKIT_URL="https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_arm64.tar.gz"; \
    else \
        SEQKIT_URL="https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz"; \
    fi && \
    wget $SEQKIT_URL -O seqkit.tar.gz && \
    tar -zxvf seqkit.tar.gz && \
    mv seqkit /usr/local/bin/ && \
    rm seqkit.tar.gz

# 安装 MAFFT
WORKDIR /home/tools/
RUN wget -c https://mafft.cbrc.jp/alignment/software/mafft-7.525-linux.tgz && \
    tar -zxvf mafft-7.525-linux.tgz && \
    mv mafft-linux64 mafft && \
    cd mafft && \
    chmod a+x ./mafft.bat && \
    mv mafft.bat mafft && \
    ln -s /home/tools/mafft/mafft /usr/local/bin/mafft && \
    cd .. && \
    rm mafft-7.525-linux.tgz

# 设置环境变量
ENV PATH="/home/tools:/usr/local/bin:${PATH}"

# 验证安装
RUN which fastp && \
    which seqkit && \
    which ustacks && \
    which mafft

# 复制文件
WORKDIR /home/app
COPY requirements.txt ./
RUN pip install -r requirements.txt
COPY ./*.py ./
COPY classification_proposal/classification.py ./
COPY *_tmp.conf ./
COPY . /home/app

# 设置工作目录
WORKDIR /home/app

RUN apt-get update && apt-get install -y \
    python3-dev \
    libfreetype6-dev \
    libpng-dev \
    pkg-config