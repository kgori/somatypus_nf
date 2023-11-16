FROM ubuntu:latest

# Install packages
RUN apt-get update && apt-get install -y --no-install-recommends build-essential python2.7-dev libbz2-dev zlib1g-dev liblzma-dev git wget curl python-pip autoconf ncurses-dev automake pkg-config && rm -rf /var/lib/apt/lists/*
RUN ln -s /usr/bin/python2.7 /usr/bin/python && python -m pip install "cython==0.29.15" && python -m pip install "pysam==0.20.0" && python -m pip install "enum34==1.1.10"

RUN git clone --recurse-submodules https://github.com/samtools/htslib.git
WORKDIR /htslib
RUN autoheader && autoreconf --install && ./configure && make && make install
WORKDIR /

RUN git clone --recurse-submodules https://github.com/samtools/samtools.git
WORKDIR /samtools
RUN autoheader && autoconf -Wno-syntax && ./configure && make && make install
WORKDIR /

RUN git clone --recurse-submodules https://github.com/samtools/bcftools.git
WORKDIR /bcftools
RUN make && make install
WORKDIR /

RUN git clone --recurse-submodules https://github.com/vcftools/vcftools.git
WORKDIR /vcftools
RUN ./autogen.sh && ./configure && make && make install
WORKDIR /

RUN rm -rf htslib && rm -rf samtools && rm -rf bcftools && rm -rf vcftools

RUN git clone --recurse-submodules https://github.com/kgori/Platypus.git
WORKDIR /Platypus
RUN git checkout container_compatible && make && echo '#!/bin/bash' > /usr/local/bin/platypus && echo "python /Platypus/bin/Platypus.py \$@" >> /usr/local/bin/platypus && chmod +x /usr/local/bin/platypus
WORKDIR /

RUN ldconfig && useradd --create-home myuser && chown -R myuser:myuser /home/myuser
USER myuser
WORKDIR /home/myuser
CMD ["bash"]
