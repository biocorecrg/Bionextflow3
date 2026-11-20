#!/usr/local/bin/Rscript --vanilla

source("/functions.R")

# Override create_pca_data so all colData(vsd) annotation columns are exported
create_pca_data <- function(vsd, condition, pcnum, colors) {

	p <- plotPCA_plus(vsd, intgroup=condition, returnData = TRUE, pcnum=pcnum)
	pca_data <- p[["data"]]
	pca_variance <- p[["variance"]]

	# Generate distinct colors dynamically for the factor levels of the condition column
	cond_vals <- as.character(pca_data[[condition]])
	cond_vals[is.na(cond_vals)] <- "NA"
	unique_vals <- unique(cond_vals)
	num_unique <- length(unique_vals)

	# Generate distinct colors dynamically (recycle if more levels than palette)
	color_palette <- rep_len(colors, num_unique)

	# Create a lookup table
	color_map <- setNames(color_palette, unique_vals)

	# Add a new column 'color' to pca_data
	pca_data$color <- unname(color_map[cond_vals])

	# Include all additional sample annotation columns from colData(vsd)
	coldata_df <- as.data.frame(colData(vsd))
	matched_idx <- match(as.character(pca_data$name), rownames(coldata_df))
	reserved_cols <- c("file", "sample", "fastq_1", "fastq_2", "file_1", "file_2", "fastq", "file_path", "read_1", "read_2", "sample_name", "sizeFactor", "replaceable", "name", "group", "color", colnames(pca_data))
	candidate_cols <- setdiff(colnames(coldata_df), reserved_cols)
	candidate_cols <- candidate_cols[!grepl("^PC[0-9]+$", candidate_cols)]

	sample_names_norm <- trimws(as.character(pca_data$name))
	extra_kept <- character(0)
	for (col_name in candidate_cols) {
		col_vals <- coldata_df[[col_name]][matched_idx]
		col_vals_norm <- trimws(as.character(col_vals))
		non_na <- !is.na(col_vals_norm) & !is.na(sample_names_norm)
		if (any(non_na) && all(col_vals_norm[non_na] == sample_names_norm[non_na])) {
			next
		}
		pca_data[[col_name]] <- col_vals
		extra_kept <- c(extra_kept, col_name)
	}

	pc_cols <- grep("^PC[0-9]+$", colnames(pca_data), value = TRUE)
	ordered_cols <- unique(c("group", condition, extra_kept, "name", pc_cols, "color"))
	ordered_cols <- c(ordered_cols, setdiff(colnames(pca_data), ordered_cols))
	pca_data <- pca_data[, ordered_cols, drop = FALSE]

	return(list(data = pca_data, variance = pca_variance))
}

# Export selected gene expression values (long TSV) for interactive MultiQC visualization
export_genes_expression <- function(expr_mat, dds, genes_desc, genes_raw, out_file, colors, condition) {
	genes_req <- trimws(as.character(genes_raw))
	genes_req <- genes_req[nzchar(genes_req)]
	if (length(genes_req) == 0) {
		return(invisible(NULL))
	}

	g_id_raw <- as.character(genes_desc$gene.id)
	g_name_raw <- as.character(genes_desc$gene.name)
	gene_name_norm <- normalize_gene_name(g_name_raw)
	gene_id_norm <- normalize_gene_name(g_id_raw)
	gene_id_novers_norm <- normalize_gene_name(sub("\\..*$", "", g_id_raw))
	expr_rownames <- rownames(expr_mat)

	matched_ids <- character(0)
	matched_names <- character(0)
	missing_queries <- character(0)

	for (q in genes_req) {
		q_norm <- normalize_gene_name(q)
		if (!nzchar(q_norm)) {
			next
		}
		idx <- which(gene_name_norm == q_norm & (g_id_raw %in% expr_rownames))
		if (length(idx) == 0) {
			idx <- which((gene_id_norm == q_norm | gene_id_novers_norm == q_norm) & (g_id_raw %in% expr_rownames))
		}
		if (length(idx) == 0) {
			missing_queries <- c(missing_queries, q)
			next
		}
		for (i in idx) {
			gid <- g_id_raw[i]
			if (!(gid %in% matched_ids)) {
				gname <- trimws(g_name_raw[i])
				if (is.na(gname) || !nzchar(gname)) {
					gname <- gid
				}
				matched_ids <- c(matched_ids, gid)
				matched_names <- c(matched_names, gname)
			}
		}
	}

	if (length(matched_ids) == 0) {
		warning(paste0("Gene names: [", paste(genes_req, collapse = ", "), "] NOT found or filtered out in ", out_file))
		return(invisible(NULL))
	}
	message(paste0("Found genes for ", out_file, ": [", paste(matched_names, collapse = ", "), "]"))
	if (length(missing_queries) > 0) {
		warning(paste0("Requested genes not found or filtered out by min_count: [", paste(missing_queries, collapse = ", "), "]"))
	}

	sub_mat <- expr_mat[matched_ids, , drop = FALSE]
	dup_names <- unique(matched_names[duplicated(matched_names)])
	display_names <- ifelse(matched_names %in% dup_names, paste0(matched_names, " (", matched_ids, ")"), matched_names)

	sample_ids <- colnames(sub_mat)
	coldata_df <- as.data.frame(colData(dds))
	matched_sidx <- match(sample_ids, rownames(coldata_df))

	if (condition %in% colnames(coldata_df)) {
		cond_vals <- as.character(coldata_df[[condition]][matched_sidx])
	} else {
		cond_vals <- rep("NA", length(sample_ids))
	}
	cond_vals[is.na(cond_vals) | !nzchar(trimws(cond_vals))] <- "NA"

	unique_vals <- unique(cond_vals)
	color_palette <- rep_len(colors, length(unique_vals))
	color_map <- setNames(color_palette, unique_vals)
	sample_colors <- unname(color_map[cond_vals])

	reserved_cols <- c("gene_id", "gene_name", "sample", "group", "value", "color", "file", "fastq_1", "fastq_2", "file_1", "file_2", "fastq", "file_path", "read_1", "read_2", "sample_name", "sizeFactor", "replaceable")
	candidate_cols <- setdiff(colnames(coldata_df), reserved_cols)
	candidate_cols <- candidate_cols[!grepl("^PC[0-9]+$", candidate_cols)]
	candidate_cols <- unique(c(intersect(condition, candidate_cols), candidate_cols))

	sample_names_norm <- trimws(as.character(sample_ids))
	annot_list <- list()
	for (col_name in candidate_cols) {
		col_vals <- coldata_df[[col_name]][matched_sidx]
		col_vals_norm <- trimws(as.character(col_vals))
		non_na <- !is.na(col_vals_norm) & !is.na(sample_names_norm)
		if (any(non_na) && all(col_vals_norm[non_na] == sample_names_norm[non_na])) {
			next
		}
		clean_vals <- ifelse(is.na(col_vals) | !nzchar(col_vals_norm), "NA", gsub("[\t\r\n]+", " ", col_vals_norm))
		annot_list[[col_name]] <- clean_vals
	}

	n_g <- length(matched_ids)
	n_s <- length(sample_ids)
	vals <- as.numeric(t(sub_mat))
	vals <- ifelse(is.finite(vals), round(vals, 4), NA_real_)

	out_df <- data.frame(
		gene_id = rep(matched_ids, each = n_s),
		gene_name = rep(display_names, each = n_s),
		sample = rep(sample_ids, times = n_g),
		group = rep(cond_vals, times = n_g),
		value = vals,
		color = rep(sample_colors, times = n_g),
		stringsAsFactors = FALSE
	)
	for (col_name in names(annot_list)) {
		out_df[[col_name]] <- rep(annot_list[[col_name]], times = n_g)
	}

	write.table(out_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
	invisible(out_df)
}

create_genes_expression_tables <- function(dds, gene_list, groups, desc_genes, prefix = "", colors, condition, vsd = NULL) {
	genes_raw <- unlist(strsplit(gene_list, ","))
	if (is.null(vsd)) {
		vsd <- makeVST(dds, FALSE)
	}
	vst_mat <- assay(vsd)

	norm <- counts(dds, normalized = TRUE)
	norm_log <- log2(norm + 1)

	pfx <- if (nzchar(prefix)) paste0(prefix, "_") else ""
	export_genes_expression(
		expr_mat = vst_mat,
		dds = dds,
		genes_desc = desc_genes,
		genes_raw = genes_raw,
		out_file = paste0(pfx, "vst_genes_expression.tsv"),
		colors = colors,
		condition = condition
	)
	export_genes_expression(
		expr_mat = norm_log,
		dds = dds,
		genes_desc = desc_genes,
		genes_raw = genes_raw,
		out_file = paste0(pfx, "log2_deseq_norm_genes_expression.tsv"),
		colors = colors,
		condition = condition
	)
}



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
parser$add_argument("-genes", "--Gene_list", type="character", default="", help="Comma separated list of genes for interactive gene expression plots between conditions (e.g: -genes Slc26a4,Cd276,Cyth4)")
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
condition<-trimws(gsub("~","",args$Assay_field))
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


##################################  Gene list expression tables for vst and log2(deseq-normalized + 1) ##################################

# Extended color palette for gene expression plots
boxplot_colors <- c("black","#1C9BCD","magenta","#E69F00","cyan","red","#14db1e","#dad60e","blue","darkred","#09e0d5","darkgreen","#bf3dbf","#8b836e","#4f0e29")

if (args$Desc_genes != "" & args$Gene_list != "") {

	create_genes_expression_tables(dds, args$Gene_list, groups, desc, prefix = "", colors = boxplot_colors, condition = condition, vsd = vsd)

}

################################## BATCH effect removal ##################################

if (args$add_batch & ("batch" %in% colnames(colData(dds)))) {
	print("Using batch column to remove batch effect with Combat-seq")
	if (grepl("batch", args$Assay_field)) {
		print("ERROR DO NOT PUT THE BATCH AS CONTROLLING FACTOR!!!")}

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


	####### Gene list expression tables for vst and log2(deseq-normalized + 1)

	if (args$Desc_genes != "" & args$Gene_list != "") {

		create_genes_expression_tables(se2, args$Gene_list, groups, desc, prefix = "batch", colors = boxplot_colors, condition = condition, vsd = vsd)
	}



}

################################## Creating GSEA files ##################################

if (args$GSEA_files) {

	makeGseaFiles(norm_counts, args$Desc_exp, condition)

}
