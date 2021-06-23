FROM nfcore/base:1.14
LABEL authors="charlotte.west@crick.ac.uk" \
        description="Docker image containing all requirements for development of thresholding peaks and plotting"

# Install conda packages
COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/reporting/bin:$PATH
