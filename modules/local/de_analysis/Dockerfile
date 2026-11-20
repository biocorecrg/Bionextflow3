# This builds ghcr.io/biocorecrg/biocorecrg_deseq
FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.authors="Luca Cozzuto <lucacozzuto@gmail.com>"


RUN Rscript -e "BiocManager::install(c('GenomicRanges', 'SummarizedExperiment', 'genefilter', 'geneplotter', 'DESeq2', 'tximport', 'pheatmap', 'GEOquery', 'Biobase', 'affy', 'limma', 'ggplot2', 'ggdendro', 'ggrepel', 'AnnotationDbi', 'hgu133plus2.db', 'hexbin', 'dendextend', 'sva', 'EnhancedVolcano'), ask=FALSE, update=FALSE)"

RUN Rscript -e "install.packages(c('openxlsx', 'argparse', 'tidyr', 'tidyverse', 'yaml'), Ncpus=4)"

RUN Rscript -e "pkgs <- c('argparse', 'stringr', 'DESeq2', 'tximport', 'pheatmap', 'yaml', 'tidyverse', 'dplyr', 'tidyr', 'ggrepel', 'sva', 'openxlsx', 'EnhancedVolcano'); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop('Missing R packages: ', paste(missing, collapse = ', ')); stopifnot(exists('ComBat_seq', where = asNamespace('sva'))); source('/functions.R'); cat('All required R packages and /functions.R verified.\n')"

COPY modules/local/de_analysis/resources/usr/bin/functions.R /functions.R
