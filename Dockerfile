# Use a specific version of R: R 4.3
FROM --platform=linux/amd64 rocker/verse:4.3.2

# 1. Install System Dependencies
# SeuratDisk requires HDF5 libraries to compile
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libglpk-dev \
    patch \
    && rm -rf /var/lib/apt/lists/*

# 2. Install Bioconductor & edgeR
# We install BiocManager first, then use it to pull edgeR
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')" \
    && R -e "BiocManager::install('edgeR', update = FALSE, ask = FALSE)"

RUN R -e "install.packages('hdf5r', configure.args = '--with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial')"

# 3. Install R Dependencies
RUN R -e "install.packages('remotes')" \
    && R -e "remotes::install_version('SeuratObject', version = '5.0.0', upgrade = 'never')" \
    && R -e "remotes::install_github('mojaveazure/seurat-disk', upgrade = 'never')"

# 4. Install package dependencies
RUN R -e "install.packages(c( 'circular', 'enrichR', 'ggpubr', 'dplyr', 'foreach','graphics','rlang','stats','matrixStats','tibble','tidyr','Matrix','MetaCycle','edgeR','ggplot2','patchwork','methods'))" 

# 5. Install Orpheus package from the local context
# This copies  current folder into the container and installs it
COPY . /home/rstudio/Orpheus
## RUN R -e "devtools::install('/home/rstudio/Orpheus', dependencies=TRUE)"