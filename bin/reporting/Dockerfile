FROM continuumio/miniconda3:4.9.2
LABEL authors="chris.cheshire@crick.ac.uk" \
      description="Docker image containing all requirements plotly dash"

# Init
RUN mkdir /app

# Install conda and pip packages
COPY ./environment.yml /
COPY ./requirements.txt /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/dash/bin:$PATH

# copy files over
COPY ./app.py /app/app.py
COPY ./lib /app/lib

ENV DASH_DEBUG_MODE True
WORKDIR /app
EXPOSE 8050
CMD ["python", "app.py"]