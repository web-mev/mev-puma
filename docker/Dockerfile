FROM debian:bullseye

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    git \
    python3 \
    python3-pip \
    unzip 

# Install all software under /opt/software:
RUN mkdir -p /opt/software 
RUN mkdir -p /opt/software/resources

# pre-install pandas to avoid issues with API incompatibility with pandas versions >2.0
RUN pip3 install pandas==1.5.3

# Install the NetZoo python library
RUN wget https://github.com/netZoo/netZooPy/archive/refs/tags/0.9.14.tar.gz \
    -P /opt/software && \
    tar xzvf /opt/software/0.9.14.tar.gz -C /opt/software && \
    pip3 install -e /opt/software/netZooPy-0.9.14/


# Copy the bucket-based file into this Docker:
# Note that for the WDL run, these files will be addressed directly.
# However, failure of this step in the Docker build (e.g. as part of CI tools)
# will indicate that the files are missing and something is wrong.
# So this serves as a pre-check.
RUN wget https://webmev-public.s3.us-east-2.amazonaws.com/puma_mirna_priors.ensg.tsv -O /opt/software/resources/puma_mirna_priors.ensg.tsv
RUN wget https://webmev-public.s3.us-east-2.amazonaws.com/puma_mirna_priors.symbol.tsv -O /opt/software/resources/puma_mirna_priors.symbol.tsv

# Add local scripts/files:
ADD puma.py /opt/software/ 
