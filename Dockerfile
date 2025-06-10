FROM rocker/r-base:4.5.0

# Debian testing provides Python 3.12 (and distutils has moved to its own wheel)
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.12 python3.12-dev python3.12-venv \
    build-essential wget curl git ca-certificates gnupg \
    libxml2-dev libssl-dev libcurl4-openssl-dev libgit2-dev \
    libpng-dev libjpeg-dev libtiff5-dev zlib1g-dev pigz less vim \
    && ln -s /usr/bin/python3.12 /usr/local/bin/python \
    && rm -rf /var/lib/apt/lists/*

# Bioconductor that matches the R in this image
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org', ask=FALSE)" \
    && Rscript -e "BiocManager::install(version='3.21', ask=FALSE, update=FALSE)" \
    && Rscript -e "BiocManager::install(c('edgeR','limma','tximport','ComplexHeatmap'), ask=FALSE, update=FALSE)"

###############################################################################
#  2.  Python 3.11 + uv (system-site install via --system)
###############################################################################

ENV CARGO_HOME=/root/.cargo
# Install uv and make it available in PATH
RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    ln -s /root/.local/bin/uv /usr/local/bin/uv




# ---- dependency layer (better cache) ----
WORKDIR /workspace
COPY ./pyproject.toml uv.lock* ./

# install all Python deps straight into the image-wide interpreter
RUN uv sync --no-cache-dir   # reproducible, lock-file aware

# project code copied afterwards so dependency cache stays hot
COPY ./main_workflow /workspace/main_workflow

###############################################################################
#  3.  Command-line genomics tools
###############################################################################
ARG KALLISTO_VER=0.51.1
RUN wget -q https://github.com/pachterlab/kallisto/releases/download/v${KALLISTO_VER}/kallisto_linux-v${KALLISTO_VER}.tar.gz && \
    tar -xzf kallisto_linux-v${KALLISTO_VER}.tar.gz && \
    mv kallisto/kallisto /usr/local/bin/ && \
    rm -rf kallisto*.tar.gz kallisto

# NCBI Entrez Direct (EDirect) – official one-liner
RUN sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)" && \
    ln -s /root/edirect*/edirect.pl /usr/local/bin/edirect

###############################################################################
#  SRA Toolkit 3.2.1  (adds fasterq-dump, prefetch, etc.)
###############################################################################
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

###############################################################################
#  4.  Final image metadata
###############################################################################
WORKDIR /workspace
CMD ["/bin/bash"]

###############################################################################
#  Build:
#     docker build -t uorca:0.2 .
#  Run:
#     docker run --rm -it -v $(pwd):/workspace uorca:0.2
###############################################################################
