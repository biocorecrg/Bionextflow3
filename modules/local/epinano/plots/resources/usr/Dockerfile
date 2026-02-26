FROM r-base:4.3.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('plyr', 'ggplot2', 'ggrepel', 'MASS', 'reshape2'), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /data

CMD ["R"]
