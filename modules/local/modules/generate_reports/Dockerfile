FROM nfcore/base:1.14
LABEL authors="chris.cheshire@crick.ac.uk" \
        description="Docker image containing all requirements for development of static reports"

# Install conda packages
COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/reporting/bin:$PATH
