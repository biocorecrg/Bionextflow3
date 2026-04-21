#!/usr/local/bin/Rscript --vanilla

source(Sys.which("functions.R"))



# Load common pars
parser <- getCommonPars()

#Define desired outputs:
#SPECIFIC FEATURES:
parser$add_argument("-input", "--Input_f", required="true", type="character", help="Path or matrix file to read counts")
parser$add_argument("-type", "--Input_type", type="character", required="true", help="Type of input: counts, salmon, star, matrix")
parser$add_argument("-label", "--Add_labels", type="character", default="", help="Column to be used for labels")
parser$add_argument("-dtrans", "--Desc_Tx", type="character", default="tx2gene.csv", help="Path to read tx2gene file (ONLY FOR SALMON) [default %(default)s]")
parser$add_argument("-strand", "--Strand_id", type="integer", default=4, help="Strand to analyze, 2: unstranded, 3: forward, 4 reverse (ONLY FOR STAR) [default %(default)s]")
parser$add_argument("-batch", "--add_batch", action='store_true', help="Remove batch effect with Combat-seq using the column batch in desc file")
parser$add_argument("-pcnum", "--Number_principal_components", type="integer", default=4, help="Number of principal components (one for each column) data to extract to the data table, [default %(default)s]")
parser$add_argument("-genes", "--Gene_list", type="character", default="", help="Comma separated list of genes (max 10), for boxplot gene expression between conditions (e.g: -genes Slc26a4,Cd276,Cyth4)")
parser$add_argument("-color", "--PCA_color",type="character", default="", help="Column in desc.txt to be used for coloring PCA plot")
parser$add_argument("-gsea", "--GSEA_files", action="store_true", help="Creates cls (samples grouped by condition column) and normalized counts table for GSEA analysis")



#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()


if (args$Input_type == "counts") {
	out <- makeDDSFromCounts(args$Input_f, args$Desc_exp, args$Assay_field)
} else if (args$Input_type == "salmon") {
	out <- makeDDSFromSalmon(args$Input_f, args$Desc_exp, args$Assay_field, args$Desc_Tx)
} else if (args$Input_type == "star") {
	out <- makeDDSFromStar(args$Input_f, args$Desc_exp, args$Assay_field, args$Strand_id)
} else if (args$Input_type == "matrix") {
	out <- makeDDSFromMatrix(args$Input_f, args$Desc_exp, args$Assay_field)
} else {
	print("please define one correct input type!")
    stop()
}


dds <- filterDDS(out$dds, args$min_count)

vsd<-makeVST(dds, FALSE)

# Removing ~ for the condition field
condition<-gsub("~","",args$Assay_field) 
print(args$Number_principal_components)

##################################  PCA #############################

pca_colors <- c("black","#1C9BCD","magenta","#E69F00","cyan","red","#14db1e","#dad60e","blue","darkred","#09e0d5","darkgreen","#bf3dbf","#8b836e","#4f0e29")

no_batch_pca <- create_pca_data(vsd, condition, args$Number_principal_components, pca_colors)

# Save PCA data and variance tables, this tables will be read by multiqc_pca plugin 
write.table(no_batch_pca$data, file=paste0("PCA_data.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(no_batch_pca$variance, file=paste0("PCA_variance.tsv"), sep="\t", quote=FALSE, row.names=FALSE)


### Create PCA colored by another column name
if (args$PCA_color != "") {
	color_condition <- args$PCA_color
	prefix <- paste0("PCA_colored_by_", args$PCA_color)
	no_batch_pca_color <- create_pca_data(vsd, color_condition, args$Number_principal_components, pca_colors)
	write.table(no_batch_pca_color$data, file=paste0(prefix, "_PCA_data.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}


if (args$Desc_genes != "") {
    print("found desc gene file")

    # Get results
    desc <- makeDesc(args$Desc_genes)
    norm_counts <- printCounts(dds, desc)
	groups <- colData(dds)[[condition]]
    
}


##################################  Sample CLustering table ################################## 

sample_clustering_matrix <- sample_clustering(vsd)
write.table(sample_clustering_matrix, file="Sample_Clustering_matrix.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


##################################  Gene list expression boxplots for vst and log2(deseq-normalized) ##################################

# Extended color palette for boxplots to support up to 50 groups
boxplot_colors <- c("black","#1C9BCD","magenta","#E69F00","cyan","red","#14db1e","#dad60e","blue","darkred","#09e0d5","darkgreen","#bf3dbf","#8b836e","#4f0e29")

if (args$Desc_genes != "" & args$Gene_list != "") {

create_genes_boxplots(dds, args$Gene_list, groups, desc, prefix = "", colors = boxplot_colors, condition = condition)

}

################################## BATCH effect removal ##################################

if (args$add_batch != "") {
    if (grepl("batch", args$Assay_field)) {
        stop("ERROR DO NOT PUT THE BATCH AS CONTROLLING FACTOR!!!")
        
    }
	# if (args$add_batch == "limma"){
	# 	mat<-assay(vsd)
	# 	mm <- model.matrix(as.formula(args$Assay_field), colData(vsd))
	# 	mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
	# 	assay(vsd) <- mat

	# 	create_pca_data(vsd, condition, args$Number_principal_components, "batch_PCA", pca_colors)

	# } 
	if (args$add_batch) {
		# Raw counts
		mat <- counts(dds)

		# Batch variable
		batch <- colData(dds)$batch

		# Biological condition to preserve
		#group <- colData(dds)$condition

		# Run ComBat_seq
		mat_corrected <- ComBat_seq(
		counts = mat,
		batch = batch)

		dds <- DESeqDataSetFromMatrix(countData = mat_corrected, colData = colData(dds), design =as.formula(args$Assay_field))

		# Recalculate VST
		se2 <- DESeq(dds)
		vsd <- makeVST(se2, FALSE)

		##### PCA for the counts with batch effect removal
		
		batch_pca <- create_pca_data(vsd, condition, args$Number_principal_components, pca_colors)
		write.table(batch_pca$data, file=paste0("batch_PCA_data.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
		write.table(batch_pca$variance, file=paste0("batch_PCA_variance.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

		##### Sample CLustering table

		sample_clustering_matrix <- sample_clustering(vsd)
		write.table(sample_clustering_matrix, file="batch_Sample_Clustering_matrix.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


		####### Gene list expression boxplots for vst and log2(deseq-normalized)

		if (args$Desc_genes != "" & args$Gene_list != "") {

			create_genes_boxplots(dds, args$Gene_list, groups, desc, prefix = "batch", colors = boxplot_colors, condition = condition)
		}



	}


}

################################## Creating GSEA files ##################################

if (args$GSEA_files) {

	makeGseaFiles(norm_counts, args$Desc_exp, condition)

}