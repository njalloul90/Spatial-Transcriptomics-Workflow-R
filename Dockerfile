# Use rocker image (R + tidyverse + system libs)
FROM rocker/tidyverse:4.4.1

LABEL maintainer="your_name <your_email>"

# Install system dependencies for spatial transcriptomics packages
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c( \
  'here', 'Matrix', 'cowplot', 'patchwork', 'jsonlite', 'scales'))"

RUN R -e "if(!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')"

RUN R -e "BiocManager::install(c( \
  'Seurat', 'SeuratDisk', 'SingleCellExperiment', 'SingleR', 'celldex', \
  'edgeR', 'limma', 'SpatialExperiment' \
), ask=FALSE, update=TRUE)"

# Copy repo into image
WORKDIR /app
COPY . /app

# Default command: run R
CMD ["R"]
