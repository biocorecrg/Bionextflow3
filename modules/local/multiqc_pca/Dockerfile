# This will do biocorecrg/deseq2:2023_1
FROM bioconductor/bioconductor_docker:RELEASE_3_14

# File Author / Maintainer
MAINTAINER Luca Cozzuto <lucacozzuto@gmail.com> 

RUN Rscript -e "BiocManager::install(c('GenomicRanges', 'SummarizedExperiment', 'genefilter', 'geneplotter'))"
RUN Rscript -e "BiocManager::install(c('DESeq2', 'tximport', 'pheatmap'), dependencies=TRUE)"
RUN Rscript -e "BiocManager::install(c('GEOquery', 'Biobase', 'affy', 'limma', 'ggplot2', 'ggdendro', 'AnnotationDbi', 'hgu133plus2.db', 'hexbin', 'dendextend','sva'), dependencies=TRUE)"

RUN Rscript -e "devtools::install_github('zhangyuqing/sva-devel')"

RUN Rscript -e "install.packages('openxlsx')"

RUN Rscript -e "BiocManager::install('EnhancedVolcano')"

RUN Rscript -e "install.packages('argparse')" 

RUN Rscript -e "install.packages('tidyr')"

RUN Rscript -e "install.packages('tidyverse')"

RUN Rscript -e "install.packages('yaml')"



