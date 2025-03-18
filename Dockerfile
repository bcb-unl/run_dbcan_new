FROM continuumio/miniconda3

# set work dir
WORKDIR /app

# install conda packages
RUN conda create -y -n run_dbcan python=3.8 dbcan -c conda-forge -c bioconda && conda clean -ya

# add run_dbcan environment to path
RUN echo "source activate run_dbcan" > ~/.bashrc
ENV PATH /opt/conda/envs/run_dbcan/bin:$PATH

# download and make the databases
RUN run_dbcan database --db_dir db


ENTRYPOINT ["run_dbcan"]

CMD [ "run_dbcan -h" ]
