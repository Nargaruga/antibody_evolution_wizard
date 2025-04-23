FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y wget

# Set up working directory
COPY . /efficient-evolution
WORKDIR /efficient-evolution

# Set up the environment
RUN conda env create -f environment.yml
RUN conda run -n efficient-evolution pip install .