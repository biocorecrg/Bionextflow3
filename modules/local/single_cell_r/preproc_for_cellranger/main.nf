process PREPROC_CELLRANGER {
    tag "$meta.id"
    label 'process_low'

    container 'biocorecrg/sc_benchmark:0.2'

    input:
    tuple val(meta), path(quants_folder)

    output:
    tuple val(meta), path("*.rds"), emit: rds

    when:
    task.ext.when == null || task.ext.when

    script:
 
    args = task.ext.args ?: ''
 
    def prefix = task.ext.prefix ?: "${meta.id}"


       """
cat > CMD.R << 'EOL'

library(dplyr)
library(Seurat)
library(patchwork)
library(tximport)
library("ggplot2")
library("SingleR")
library("celldex")
library("stringr")
 
data_dir <- "${quants_folder}/filtered_feature_bc_matrix"
data <- Read10X(data.dir = data_dir)


if(is.list(data)) {
	olddata <- data
	data <- olddata[["Gene Expression"]]
}

seurObj = CreateSeuratObject(counts = data)
# save seurat object
saveRDS(seurObj, file = "${prefix}.rds")
quit("no")

EOL

	Rscript CMD.R 
        """
}
