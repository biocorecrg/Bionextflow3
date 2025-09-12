process PREPROC_PARSE {
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

expression_matrix <- ReadParseBio("${quants_folder}/DGE_filtered")
gene_genomes <- read.csv(paste0("${quants_folder}/DGE_filtered", "/all_genes.csv"))

# For multi species... 
#comb_gene<-paste(gene_genomes\$genome, rownames(expression_matrix), sep="-")

# if empty gene names are present, name them unknown.
rownames(expression_matrix)[rownames(expression_matrix) == ""] <- "unknown"
# Read in cell meta data
cell_meta <- read.csv(paste0("${quants_folder}/DGE_filtered", "/cell_metadata.csv"), row.names = 1)
seurObj <- CreateSeuratObject(expression_matrix, names.field = 0, meta.data = cell_meta, min.cells = 3, min.features = 200, project = "${prefix}")
Idents(seurObj) <- seurObj@meta.data\$orig.ident


# save seurat object
saveRDS(seurObj, file = "${prefix}.rds")
quit("no")

EOL

	Rscript CMD.R 
        """
}
