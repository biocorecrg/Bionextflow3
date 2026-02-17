/*
*  Seurat module
*/

process SEURAT {
    tag "$meta.id"
    label 'process_high_memory'

    container 'docker.io/biocorecrg/sc_benchmark:0.2'

    input:
    tuple val(meta), path("input_ori.rds")
    val(genome)
    path(rrna_genes)

    output:
    tuple val(meta), path("*.pdf"), emit: pdfs
    tuple val(meta), path("*.rds"), emit: rds
    path "versions.yml"           , emit: versions



    when:
    task.ext.when == null || task.ext.when

    script:
 
    args = task.ext.args ?: ''
 
    def prefix = task.ext.prefix ?: "${meta.id}"

    def globalsMax = (task.memory.toBytes() / 2) as long

    def script_anno = ""
    
    if (genome == "human" || genome == "mouse") {
    	if (genome == "human") {
			script_anno = "ref.data <- HumanPrimaryCellAtlasData()"
    	} else if (genome == "mouse") {
        	script_anno = "ref.data <- MouseRNAseqData()"
    	}
    	script_anno = script_anno + """    
sce <- as.SingleCellExperiment(DietSeurat(seurObj))

predictions.main <- SingleR(test=sce, assay.type.test=1, 
    ref=ref.data, labels=ref.data\$label.main)

table(predictions.main\$labels)

seurObj@meta.data\$predictions.main <- predictions.main\$labels

seurObj <- SetIdent(seurObj, value = "predictions.main")

pdf(paste("${prefix}", "_ann.pdf", sep=""), width=10, height=10)
DimPlot(seurObj, label = T , repel = T, label.size = 3) + NoLegend()
dev.off()
"""
	} 	

    def vinplot_cmd  = """VlnPlot(seurObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)"""

   if (rrna_genes) {
    vinplot_cmd = """
    rRNA_genes <- readLines(\"${rrna_genes}\")
    print(paste0("Searching for ", length(rRNA_genes), " rRNA genes"))
    # Filtra quelli presenti nel Seurat object
    rRNA_genes<-gsub("_", "-", rRNA_genes)
    rRNA_genes_filtered <- rRNA_genes[rRNA_genes %in% rownames(seurObj)]
    print(paste0("Found ", length(rRNA_genes_filtered), " rRNA genes in the seurat object"))
    seurObj[["percent.rRNA"]] <- PercentageFeatureSet(seurObj, features = rRNA_genes_filtered)
	VlnPlot(seurObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rRNA"), ncol = 4, pt.size = 0)
"""
   }

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
 
options(future.globals.maxSize = ${globalsMax} ) 

seurObj <-readRDS("input_ori.rds")

# Extract mitochondrial genes (case-insensitive)
mt_genes <- grep("^mt-", rownames(seurObj[["RNA"]]), value = TRUE, ignore.case = TRUE)

# Check if any mitochondrial genes were found
if (length(mt_genes) == 0) {
  stop("No mitochondrial genes (starting with 'mt-' or 'MT-'') found in the dataset.")
}

# Compute percent mitochondrial content
seurObj[["percent.mt"]] <- PercentageFeatureSet(seurObj, features = mt_genes, assay = 'RNA')

pdf(paste("${prefix}", "_vp.pdf", sep=""), width=10)
${vinplot_cmd}
dev.off()

# subsetting and normalize 
# We subset cells with less than 0.05 percentile MT and between 200 and 0.99 percentile of nFeatures 
cutoff.mt<-round(quantile(seurObj[["percent.mt"]][, 1], c(.95)))[[1]]
cutoff.nFeat<-round(quantile(seurObj[["nFeature_RNA"]][, 1], c(.99)))[[1]]



pdf(paste("${prefix}", "_fc.pdf", sep=""), width=10)
plot1 <- FeatureScatter(seurObj, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position="none") + geom_hline(yintercept = cutoff.mt) + annotate("text", x=-100, y=cutoff.mt+2, label= cutoff.mt)
plot2 <- FeatureScatter(seurObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position="none") + geom_hline(yintercept = 200) + geom_hline(yintercept = cutoff.nFeat) + annotate("text", x=c(-400,-400), y=c(300, cutoff.nFeat+100), label= c(200, cutoff.nFeat))
plot1 + plot2
dev.off()

print(paste("Number of cells before subsetting:", ncol(seurObj)))

print(paste0("Subsetting cells using the following criteria: nFeature_RNA > 200 & nFeature_RNA < ", cutoff.nFeat[1], " & percent.mt < ",cutoff.mt))
seurObj <- subset(seurObj, subset = nFeature_RNA > 200 & nFeature_RNA < cutoff.nFeat[1] & percent.mt < cutoff.mt)

print(paste("Number of cells after subsetting:", ncol(seurObj)))

seurObj <- SCTransform(seurObj)
seurObj <- FindVariableFeatures(seurObj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurObj), 10)

pdf(paste("${prefix}", "_vf.pdf", sep=""), width=10)
plot1 <- VariableFeaturePlot(seurObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# PCA
seurObj <- RunPCA(seurObj, features = VariableFeatures(object = seurObj))

pdf(paste("${prefix}", "_ep.pdf", sep=""), width=10)
ElbowPlot(seurObj)
dev.off()

# Here we use 1:15 to PCs. So it might be wise to have a look at elbowplot 
# before trusting these results!

seurObj <- FindNeighbors(seurObj, dims = 1:15)
seurObj <- FindClusters(seurObj, resolution = 0.5)

seurObj <- RunUMAP(seurObj, dims = 1:15)

pdf(paste("${prefix}", "_dp.pdf", sep=""), width=10)
DimPlot(seurObj, reduction = "umap")
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive
# ones

seurObj.markers <- FindAllMarkers(seurObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurObj.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(seurObj, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

seurObj.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

pdf(paste("${prefix}", "_hm.pdf", sep=""), width=10, height=20)
DoHeatmap(seurObj, features = top10\$gene) + NoLegend()
dev.off()


${script_anno}

# save seurat object
saveRDS(seurObj, file = "${prefix}.rds")

out_file <- "versions.yml"

versions <- c(
  "    ${task.process}:",
  paste("        seurat:", packageVersion("Seurat")),
  paste("        singler:", packageVersion("SingleR")),
  paste("        celldex:", packageVersion("celldex"))
)

writeLines(versions, con = out_file)
quit("no")

EOL

	Rscript CMD.R 



	"""

}

 
