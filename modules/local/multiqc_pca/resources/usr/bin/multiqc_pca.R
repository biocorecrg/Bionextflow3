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
parser$add_argument("-batch", "--add_batch", type="character", default="", help="Remove batch effect using the column batch in desc file: limma or combat")
parser$add_argument("-pcnum", "--Number_principal_components", type="integer", default=4, help="Number of principal components (one for each column) data to extract to the data table, [default %(default)s]")
parser$add_argument("-genes", "--Gene_list", type="character", default="", help="Comma separated list of genes (max 10), for boxplot gene expression between conditions (e.g: -genes Slc26a4,Cd276,Cyth4)")


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


pca_colors <- c("black","#1C9BCD","magenta","cyan","red","blue","darkred","#E69F00","#db1432","darkgreen","#14db1e","#db9c14","#3dbfac","#5c5952","#4f0e29")

# Extended color palette for boxplots to support up to 50 groups
boxplot_colors <- c("black","#1C9BCD", "darkred","darkgreen", "#E69F00", "blue", "magenta")

create_pca_data(vsd, condition, args$Number_principal_components, "PCA", pca_colors)


if (args$Desc_genes != "") {
    print("found desc gene file")

    # Get results
    desc <- makeDesc(args$Desc_genes)
    norm_counts <- printCounts(dds, desc)
	groups <- colData(dds)[[condition]]
    
}

###### Sample CLustering table

sample_clustering_matrix <- sample_clustering(vsd)
write.table(sample_clustering_matrix, file="Sample_Clustering_matrix.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


####### Gene list expression boxplots for vst and log2(deseq-normalized)

if (args$Desc_genes != "" & args$Gene_list != "") {

create_genes_boxplots(dds, args$Gene_list, groups, desc, prefix = "", colors = boxplot_colors, condition = condition)

}

### Creating files for BATCH effect removal 

if (args$add_batch != "") {
    if (grepl("batch", args$Assay_field)) {
        stop("ERROR DO NOT PUT THE BATCH AS CONTROLLING FACTOR!!!")
        
    }
	if (args$add_batch == "limma"){
		mat<-assay(vsd)
		mm <- model.matrix(as.formula(args$Assay_field), colData(vsd))
		mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
		assay(vsd) <- mat

		create_pca_data(vsd, condition, args$Number_principal_components, "batch_PCA", pca_colors)

	} else if (args$add_batch == "combat") {
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

		create_pca_data(vsd, condition, args$Number_principal_components, "batch_PCA", pca_colors)

		##### Sample CLustering table

		sample_clustering_matrix <- sample_clustering(vsd)
		write.table(sample_clustering_matrix, file="batch_Sample_Clustering_matrix.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


		####### Gene list expression boxplots for vst and log2(deseq-normalized)

		if (args$Desc_genes != "" & args$Gene_list != "") {

			create_genes_boxplots(dds, args$Gene_list, groups, desc, prefix = "batch", colors = boxplot_colors, condition = condition)
		}



	}


}


