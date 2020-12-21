FROM ubuntu:bionic
RUN apt-get update \
&& apt-get install -y libboost-all-dev \
&& apt-get install -y zlib1g-dev \
&& apt-get install -y gengetopt \
&& apt-get install -y gnulib \
&& apt-get install -y wget \
&& apt-get install -y git

RUN mkdir /data && cd /data \
&& git clone https://github.com/Lexogen-Tools/idemuxcpp \
&& cd idemuxcpp \
&& autoreconf -i && ./configure && make && make install \
&& cd /data

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
&& mkdir /root/.conda \
&& bash Miniconda3-latest-Linux-x86_64.sh -b \
&& rm -f Miniconda3-latest-Linux-x86_64.sh 

RUN cd /data \
&& git clone -b idemuxCPP_bowtie2 http://tdrozd:PCVsDSegxZau9pyzgcab@10.90.1.56:10080/Bioinfo/sars-cov-2_analysis.git \
&& cd sars-cov-2_analysis && conda env create -f environment.yml

ENV PATH /root/miniconda3/envs/lex_sars-cov-2/bin:${PATH}

RUN /bin/bash -c "source activate lex_sars-cov-2"

WORKDIR /data/sars-cov-2_analysis
