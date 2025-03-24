FROM mambaorg/micromamba:2-ubuntu22.04 AS source

USER $MAMBA_USER

# set work directory
WORKDIR /tmp

# copy environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# create environment from yml file
RUN micromamba env create -f environment.yml \
  && micromamba clean --all --yes \
  && rm /tmp/environment.yml

FROM mambaorg/micromamba:2-ubuntu22.04 AS final

ARG BUILD_VERSION=v5.0.0

################## METADATA ######################
LABEL base_image="micromamba:2-ubuntu22.04"
LABEL version="1"
LABEL software="dbcan and dbcan utils"
LABEL software.version="${BUILD_VERSION}"
LABEL about.summary="Docker image for dbcan and dbcan utils"
LABEL about.tags="annotation, metagenomics, dbcan"
LABEL about.license="MIT"

COPY --from=source /opt/conda/envs/run_dbcan /opt/conda/envs/run_dbcan

# Won't be activated unless you we login via /bin/bash command
# RUN echo "micromamba activate fastcat" >> ~/.bashrc

# set fastcat as the default environment
# !!! this step is important or you cannot packages installed in this environment
ENV ENV_NAME run_dbcan
ENV PATH=/opt/conda/envs/run_dbcan/bin:$PATH

WORKDIR /app

# default command
CMD ["run_dbcan", "--help"]

# add mantainer
MAINTAINER  Xinpeng Zhang <xzhang55@huskers.unl.edu>, Haidong Yi <hyi@stjude.org>
