FROM rocker/r-base:4.5.0

# Install system dependencies including Python 3.13
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.13 python3.13-dev python3.13-venv \
    build-essential cmake wget curl git ca-certificates gnupg \
    libxml2-dev libssl-dev libcurl4-openssl-dev libgit2-dev \
    libpng-dev libjpeg-dev libtiff5-dev zlib1g-dev pigz less vim \
    hdf5-tools libhdf5-dev hdf5-helpers libhdf5-serial-dev \
    autoconf \
    && ln -s /usr/bin/python3.13 /usr/local/bin/python \
    && rm -rf /var/lib/apt/lists/*

# Install uv
ENV CARGO_HOME=/root/.cargo
RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    ln -sf /root/.local/bin/uv /usr/local/bin/uv

# Install R packages
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org', ask=FALSE)" \
    && Rscript -e "BiocManager::install(version='3.21', ask=FALSE, update=FALSE)" \
    && Rscript -e "BiocManager::install(c('edgeR','limma','tximport','ComplexHeatmap','gplots'), ask=FALSE, update=FALSE)"

# Install genomics tools - compile Kallisto from source following official instructions
WORKDIR /tmp
RUN git clone https://github.com/pachterlab/kallisto.git && \
    cd kallisto && \
    cd ext/htslib && \
    autoheader && \
    autoconf && \
    cd ../.. && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    cd / && \
    rm -rf /tmp/kallisto

# Install SRA Toolkit
ARG SRA_VER=3.2.1
RUN set -eux; \
    arch="$(uname -m)"; \
    case "$arch" in \
    x86_64)   SRA_ARCH=ubuntu64 ;; \
    aarch64|arm64) SRA_ARCH=linux-aarch64 ;; \
    *) echo "Unsupported arch $arch"; exit 1 ;; \
    esac; \
    wget -q "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VER}/sratoolkit.${SRA_VER}-${SRA_ARCH}.tar.gz"; \
    tar -xzf "sratoolkit.${SRA_VER}-${SRA_ARCH}.tar.gz" -C /opt; \
    rm  "sratoolkit.${SRA_VER}-${SRA_ARCH}.tar.gz"; \
    ln -s /opt/sratoolkit.${SRA_VER}-${SRA_ARCH}/bin/* /usr/local/bin/; \
    mkdir -p /root/.ncbi && \
    printf '/LIBS/GUID = "docker-build-guid"\nconfig/default = "true"\n' > /root/.ncbi/user-settings.mkfg

# Install NCBI Entrez Direct (EDirect)
RUN set -eux; \
    echo ">>> Installing EDirect …"; \
    sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"; \
    find /root/edirect -maxdepth 1 -type f -perm -u+x -exec ln -sf {} /usr/local/bin/ \; ; \
    /root/edirect/esearch -version | head -1

# Set up workspace and create proper venv
WORKDIR /workspace

# Copy dependency files first for better caching
COPY pyproject.toml uv.lock ./

# Create venv with correct Python interpreter that will be available at runtime
RUN uv venv .venv --python=/usr/bin/python3.13

# Install dependencies into the venv
RUN uv sync --no-cache

# Set up environment to use the venv
ENV PATH="/workspace/.venv/bin:${PATH}"
ENV VIRTUAL_ENV="/workspace/.venv"
ENV UV_NO_SYNC=1

# Copy project code (all migrated into uorca/)
COPY uorca/ ./uorca/

# Install the UORCA package in development mode so entry points are available
RUN uv pip install -e .

# Verify installation
RUN which uorca && uorca --help

# Final environment check
RUN echo ">>> FINAL ENVIRONMENT:"; \
    echo "PATH=$PATH"; \
    echo "VIRTUAL_ENV=$VIRTUAL_ENV"; \
    which python uv uorca edirect || true; \
    python --version; \
    uv --version

CMD ["/bin/bash"]
