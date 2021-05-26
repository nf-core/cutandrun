# Set the base image to Ubuntu 16.04
FROM rocker/rstudio:4.0.3

LABEL authors="chris.cheshire@crick.ac.uk" \
        description="Docker image containing Rstudio and a conda env for the deseq2_diff.r script"

#Install conda
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

# Update conda
RUN conda update -n base -c defaults conda

# Install conda packages
COPY environment.yml /
RUN conda env update -f /environment.yml && conda clean -a
