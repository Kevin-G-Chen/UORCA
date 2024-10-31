FROM python:3.12-slim

# Set environment variables for user creation and PATH
ENV NB_USER=myuser \
    NB_UID=1001 \
    NB_GID=1001 \
    HOME=/home/myuser \
    PATH="/usr/local/edirect:/usr/local/bin:${PATH}"

USER root

# Install necessary tools and JupyterLab dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    curl \
    build-essential \
    nodejs \
    npm \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libtiff-dev \
    libjpeg-dev \
    libpng-dev \
    libgit2-dev \
    libzmq3-dev \
    libhdf5-dev \
    libglpk-dev \
    libcairo2-dev \
    libxt-dev \
    pandoc \
    git \
    cmake \
    r-base \
    wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install JupyterLab and other necessary Python packages
RUN pip install --no-cache-dir \
    jupyterlab \
    jupyter-ai \
    jupyter-ai-magics \
    langchain_openai \
    bash_kernel && \
    python -m bash_kernel.install

# Install Kallisto
RUN wget https://github.com/pachterlab/kallisto/archive/refs/tags/v0.51.0.tar.gz -O /tmp/kallisto.tar.gz && \
    tar -zxf /tmp/kallisto.tar.gz -C /opt && \
    cd /opt/kallisto-0.51.0 && \
    mkdir build && cd build && \
    cmake .. && make && \
    ln -s /opt/kallisto-0.51.0/build/src/kallisto /usr/local/bin/kallisto && \
    rm /tmp/kallisto.tar.gz

# Build and install SRA Toolkit
RUN TMPDIR=$(mktemp -d) && \
    git clone --depth 1 --branch 3.1.1 https://github.com/ncbi/ncbi-vdb.git ${TMPDIR}/ncbi-vdb && \
    git clone --depth 1 --branch 3.1.1 https://github.com/ncbi/sra-tools.git ${TMPDIR}/sra-tools && \
    mkdir -p ${TMPDIR}/build && cd ${TMPDIR}/build && \
    cmake -S ${TMPDIR}/ncbi-vdb -B ncbi-vdb && \
    cmake --build ncbi-vdb && \
    cmake -D VDB_LIBDIR="${TMPDIR}/build/ncbi-vdb/lib" \
    -D CMAKE_INSTALL_PREFIX="/opt/sratoolkit" \
    -S ${TMPDIR}/sra-tools \
    -B sra-tools && \
    cmake --build sra-tools --target install && \
    ln -s /opt/sratoolkit/bin/* /usr/local/bin/ && \
    rm -rf ${TMPDIR}

# Install Poetry via pip
RUN pip install --no-cache-dir poetry

# Install R packages and IRkernel
RUN R -e "install.packages('pak', repos='https://cloud.r-project.org/')" && \
    R -e "pak::pkg_install('IRkernel')" && \
    R -e "IRkernel::installspec(user = FALSE)"

# Create a new user with specified UID and GID
RUN groupadd -g $NB_GID $NB_USER && \
    useradd -m -s /bin/bash -N -u $NB_UID -g $NB_GID $NB_USER

# Copy project files and set ownership using --chown
COPY --chown=$NB_USER:$NB_GID . /home/$NB_USER/work

# Copy update_packages.R and set execute permissions
COPY --chown=$NB_USER:$NB_GID update_packages.R /home/$NB_USER/update_packages.R
RUN chmod +x /home/$NB_USER/update_packages.R

# Create R_libs directory
RUN mkdir -p /home/$NB_USER/R_libs

# Set the working directory
WORKDIR /home/$NB_USER/work

# Copy packages.txt and set ownership
COPY --chown=$NB_USER:$NB_GID packages.txt /home/$NB_USER/work/packages.txt

# Install Entrez Direct (EDirect) globally
RUN curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -o /tmp/install-edirect.sh && \
    sh /tmp/install-edirect.sh -y && \
    mv /home/myuser/edirect /usr/local/ && \
    ln -sf /usr/local/edirect/* /usr/local/bin/ && \
    rm /tmp/install-edirect.sh

# Install dependencies using Poetry and set up Jupyter kernel
RUN poetry install && \
    poetry run python -m ipykernel install --user --name=poetry-env --display-name "Python (Poetry)"

# Expose Jupyter Notebook port
EXPOSE 8888

# Start JupyterLab from within the Poetry environment and apply the theme
CMD ["sh", "-c", "\
    if [ -f '/home/myuser/work/packages.txt' ]; then \
    Rscript -e \"if (!requireNamespace('pak', quietly = TRUE)) install.packages('pak', repos='https://cloud.r-project.org/'); library(pak); package_names <- unique(c(readLines('/home/myuser/work/packages.txt'), 'tidyverse')); message(package_names); pak::pkg_install(package_names)\"; \
    fi; \
    trap 'Rscript /home/myuser/update_packages.R' EXIT; \
    poetry run jupyter lab --ip=0.0.0.0 --no-browser --allow-root"]
