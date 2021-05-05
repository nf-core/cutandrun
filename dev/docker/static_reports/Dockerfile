FROM nfcore/base:1.12.1
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements for generating static reports"

# Install conda packages
COPY ./environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/reporting/bin:$PATH