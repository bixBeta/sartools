FROM rocker/r-ver:4.4.2

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    imagemagick \
    openssh-client \
    && rm -rf /var/lib/apt/lists/*

# Quarto (pinned for reproducibility; detects arm64 vs amd64 automatically)
ARG QUARTO_VERSION=1.5.57
RUN ARCH=$(dpkg --print-architecture) \
    && wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${ARCH}.deb \
    && dpkg -i quarto-${QUARTO_VERSION}-linux-${ARCH}.deb \
    && rm quarto-${QUARTO_VERSION}-linux-${ARCH}.deb

# Use Posit Package Manager for pre-built Linux binaries (Ubuntu 24.04 noble)
# This avoids source compilation failures (cmake, boost, etc.)
ENV RSPM="https://packagemanager.posit.co/cran/__linux__/noble/latest"

# CRAN packages
RUN R -e "install.packages(c('BiocManager','dplyr','progress','tibble','tidyr',\
'ggplot2','ggrepel','matrixStats','DT','viridisLite','reshape2',\
'RColorBrewer','knitr','rmarkdown','amap','xtable',\
'kableExtra','GGally','ggdendro','gridExtra'), \
repos=Sys.getenv('RSPM'), Ncpus=4)"

# Bioconductor dependencies required by SARTools 1.7.4
RUN R -e "options(repos = Sys.getenv('RSPM')); \
BiocManager::install(c('DESeq2','edgeR','biomaRt',\
'genefilter','limma'), ask=FALSE, Ncpus=4)"

# Install pinned SARTools 1.7.4 from local tarball
COPY SARTools-1.7.4.tar.gz /tmp/
RUN R CMD INSTALL /tmp/SARTools-1.7.4.tar.gz \
    && rm /tmp/SARTools-1.7.4.tar.gz
