FROM python:3.12-slim

USER root

# Install necessary tools and JupyterLab dependencies
RUN apt-get update
RUN apt-get install -y curl build-essential nodejs npm libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libfreetype6-dev libfribidi-dev libharfbuzz-dev libtiff-dev libjpeg-dev libpng-dev libtiff-dev libgit2-dev libzmq3-dev libhdf5-dev libglpk-dev libcairo2-dev libxt-dev
RUN apt-get install -y pandoc git cmake
RUN apt-get install -y r-base wget
#RUN apt-get install -y texlive-xetex texlive-fonts-recommended texlive-plain-generic wget gzip
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
# Install JupyterLab and other necessary Python packages
RUN pip install --no-cache-dir jupyterlab
RUN pip install --no-cache-dir jupyter-ai
RUN pip install --no-cache-dir jupyter-ai-magics
RUN pip install --no-cache-dir langchain_openai
RUN pip install --no-cache-dir bash_kernel && python -m bash_kernel.install

RUN wget https://github.com/pachterlab/kallisto/archive/refs/tags/v0.51.0.tar.gz -O /tmp/kallisto.tar.gz && \
    tar -zxf /tmp/kallisto.tar.gz -C /opt && \
    cd /opt/kallisto-0.51.0 && \
    mkdir build && cd build && \
    cmake .. && make && \
    ln -s /opt/kallisto-0.51.0/build/src/kallisto /usr/local/bin/kallisto && \
    rm /tmp/kallisto.tar.gz


# Build and install SRA Toolkit
RUN TMPDIR=$(mktemp -d) && \
    cd ${TMPDIR} && \
    git clone https://github.com/ncbi/ncbi-vdb.git && \
    git clone https://github.com/ncbi/sra-tools.git && \
    mkdir build && cd build && \
    cmake -S "$(cd ../ncbi-vdb; pwd)" -B ncbi-vdb && \
    cmake --build ncbi-vdb && \
    cmake -D VDB_LIBDIR="${PWD}/ncbi-vdb/lib" -D CMAKE_INSTALL_PREFIX="/opt/sratoolkit" -S "$(cd ../sra-tools; pwd)" -B sra-tools && \
    cmake --build sra-tools --target install && \
    ln -s /opt/sratoolkit/bin/* /usr/local/bin/ && \
    rm -rf ${TMPDIR}

# Install Poetry via pip
RUN pip install --no-cache-dir poetry

RUN R -e "install.packages('pak')"
RUN R -e "pak::pkg_install('IRkernel')"
RUN R -e "IRkernel::installspec(user = FALSE)"

# Create a new user
ENV NB_USER=myuser \
    NB_UID=1001 \
    NB_GID=1001 \
    HOME=/home/myuser

RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER

# Copy project files
COPY . /home/$NB_USER/work

COPY update_packages.R /home/$NB_USER/update_packages.R
RUN chmod +x /home/$NB_USER/update_packages.R

RUN mkdir -p /home/$NB_USER/R_libs

# Set the working directory
WORKDIR /home/$NB_USER/work

COPY packages.txt /home/$NB_USER/work/packages.txt

# Change ownership of the directories

RUN chown -R $NB_USER:$NB_GID /home/$NB_USER && \
    mkdir -p /home/$NB_USER/.local/share/jupyter && \
    chown -R $NB_USER:$NB_GID /home/$NB_USER/.local/share/jupyter && \
    mkdir -p /home/$NB_USER/.jupyter && \
    chown -R $NB_USER:$NB_GID /home/$NB_USER/.jupyter && \
    mkdir -p /home/$NB_USER/.cache/jupyter && \
    chown -R $NB_USER:$NB_GID /home/$NB_USER/.cache/jupyter && \
    mkdir -p /home/$NB_USER/.local/share/jupyter/lab && \
    chown -R $NB_USER:$NB_GID /home/$NB_USER/.local/share/jupyter/lab && \
    mkdir -p /home/$NB_USER/.cache && \
    chown -R $NB_USER:$NB_GID /home/$NB_USER/.cache && \
    mkdir -p /home/$NB_USER/.local/share/Trash && \
    chown -R $NB_USER:$NB_GID /home/$NB_USER/.local/share/Trash

# Ensure Poetry is accessible by setting the PATH for the new user
RUN echo 'export PATH="/usr/local/bin:$PATH"' >> /home/$NB_USER/.bashrc
RUN echo 'export PATH="/usr/local/bin:$PATH"' >> /home/$NB_USER/.profile

# Switch to the new user
USER $NB_USER

# Install Entrez Direct (EDirect)
RUN sh -c "echo y | curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh" && \
    echo "export PATH=\$HOME/edirect:\$PATH" >> /home/$NB_USER/.bashrc && \
    echo "export PATH=\$HOME/edirect:\$PATH" >> /home/$NB_USER/.profile

# Install dependencies using Poetry and set up Jupyter kernel
RUN poetry install && \
    poetry run python -m ipykernel install --user --name=poetry-env --display-name "Python (Poetry)"

# Expose Jupyter Notebook port
EXPOSE 8888

# Start JupyterLab from within the Poetry environment and apply the theme
CMD ["sh", "-c", "\
    if [ -f '/home/myuser/work/packages.txt' ]; then \
    Rscript -e \"if (!requireNamespace('pak', quietly = TRUE)) install.packages('pak'); library(pak); package_names <- unique(c(readLines('/home/myuser/work/packages.txt'), 'tidyverse'));message(package_names);pak::pkg_install(package_names)\"; \
    fi; \
    trap 'Rscript /home/myuser/update_packages.R' EXIT; \
    poetry run jupyter lab --ip=0.0.0.0 --no-browser --allow-root"]
