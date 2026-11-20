# This builds ghcr.io/biocorecrg/biocorecrg_deseq
FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.authors="Luca Cozzuto <lucacozzuto@gmail.com>"

RUN Rscript -e "options(timeout = 300); BiocManager::install(c('GenomicRanges', 'SummarizedExperiment', 'annotate', 'genefilter', 'geneplotter', 'DESeq2', 'tximport', 'pheatmap', 'GEOquery', 'Biobase', 'affy', 'limma', 'ggplot2', 'ggdendro', 'ggrepel', 'AnnotationDbi', 'hgu133plus2.db', 'hexbin', 'dendextend', 'sva', 'EnhancedVolcano'), ask=FALSE, update=FALSE)"

RUN Rscript -e "options(timeout = 300); install.packages(c('openxlsx', 'argparse', 'tidyr', 'tidyverse', 'yaml'), Ncpus=4)"

COPY modules/local/de_analysis/resources/usr/bin/functions.R /functions.R

RUN Rscript -e "pkgs <- c('argparse', 'stringr', 'DESeq2', 'tximport', 'pheatmap', 'yaml', 'tidyverse', 'dplyr', 'tidyr', 'ggrepel', 'sva', 'openxlsx', 'EnhancedVolcano'); for (p in pkgs) library(p, character.only = TRUE); stopifnot(exists('ComBat_seq', where = asNamespace('sva'))); source('/functions.R'); cat('All required R packages and /functions.R verified.\n')"
