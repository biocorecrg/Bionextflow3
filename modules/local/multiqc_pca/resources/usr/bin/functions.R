suppressMessages(library("argparse"))
suppressMessages(library("stringr"))
suppressMessages(library("DESeq2"))
#suppressMessages(library("EnhancedVolcano"))
suppressMessages(library("tximport"))
suppressMessages(library("pheatmap"))
suppressMessages(library("yaml"))
suppressMessages(library("tidyverse"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("ggrepel"))




##Argument parser:
#Create parser object

# functions
#Get common pars
getCommonPars <- function(desc_exp) {
	parser <- ArgumentParser()

	#Define desired outputs:
	#GLOBAL FEATURES:
	parser$add_argument("-desc", "--Desc_exp", default="desc.txt", type="character", help="File describing experiments, [default %(default)s]")
	parser$add_argument("-dgenes", "--Desc_genes", default="gene_desc.txt", type="character", help="Files describing genes, [default %(default)s]")
    parser$add_argument("-assay", "--Assay_field", required="true", type="character", help="Formula to be used for DE analysis, example \"~ time\"")
	parser$add_argument("-min", "--min_count", type="integer", default=1, help="Number of minimum read per row [default %(default)s]")
	return (parser)
}


makeDDSFromCounts <- function(input, desc, field) {
	# Read counts and make count table
	ff <- list.files( path = input, recursive=F, pattern = "*.counts$", full.names = TRUE )
	counts.files <- lapply( ff, read.table)
	counts <- as.data.frame( sapply( counts.files, function(x) x[ , 2 ] ) )
	fn <- basename(ff)
	fn <- gsub( ".counts", "", fn)
	colnames(counts) <- fn
	row.names(counts) <- counts.files[[1]]$V1

	#Make coldata for DESEq2
	coldata <- makeColData(desc, fn)

	#DESEQ2
	dds <- DESeqDataSetFromMatrix(countData = counts,
                    colData = coldata,
                    design = as.formula(field))

	return(list("dds" = dds, "coldata" = coldata))
}

makeDDSFromSalmon <- function(input, desc, field, desctx) {
	# List the quantification files from Salmon: one quant.sf file per sample
	files <- dir(input, recursive=TRUE, pattern="quant.sf", full.names=TRUE)
	names(files) <-basename(dirname(files))

	tx2gene <- read.table(desctx, 
		sep="\t",
		header=F)

	txi <- tximport(files,
		type = "salmon", 
		tx2gene = tx2gene)

	#Make coldata for DESEq2
	coldata <- makeColData(desc, names(files))

	#DESEQ2
	dds <- DESeqDataSetFromTximport(txi,
                        colData = coldata,
                        design = as.formula(field))
	return(list("dds" = dds, "coldata" = coldata))
}

makeDDSFromStar <- function(input, desc, field, number) {

	# Read counts and make count table
	ff <- list.files( path = input, recursive=F, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
	counts.files <- lapply( ff, read.table, skip = 4 )
	counts <- as.data.frame( sapply( counts.files, function(x) x[ , number ] ) )
	fn <- basename(ff)
	fn <- gsub("\\.ReadsPerGene\\.out\\.tab", "", fn)
	colnames(counts) <- fn
	row.names(counts) <- counts.files[[1]]$V1

	#Make coldata for DESEq2
	coldata <- makeColData(desc, fn)

	#DESEQ2
	dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = as.formula(field))
	return(list("dds" = dds, "coldata" = coldata))
}

makeDDSFromMatrix <- function(input, desc, field) {

	if (basename(input) == "annotated_counts.txt.gz") {
 	 	counts <- read.delim(gzfile(input), row.names = 4, check.names = FALSE)
	} else {
 	 	counts <- read.csv(input, row.names = 1, check.names = FALSE)
    }
    
	if ('gene.name' %in% names(counts)) {
		counts$gene.name<-NULL
		counts$gene.type<-NULL
	}
	
	else if ('Strand' %in% names(counts)) {
		counts$Chr<-NULL
		counts$Start<-NULL
		counts$End<-NULL
		counts$Score<-NULL
		counts$Strand<-NULL
		counts$gene.id<-NULL
		counts$gene.name<-NULL
		counts$gene.type<-NULL
	}


    
	#Make coldata for DESEq2
	coldata <- makeColData(desc, names(counts))

	#DESEQ2
	dds <- DESeqDataSetFromMatrix(countData = counts,
                        colData = coldata,
                        design = as.formula(field))

	return(list("dds" = dds, "coldata" = coldata))
}


#Make coldata for DESEq2
makeColData <- function(desc_exp, fn) {
	desc_data<-read.table(desc_exp, sep="\t", header=T)
	coldata<-desc_data[match(fn, desc_data$file),]
	row.names(coldata)<-fn
	if ("group" %in% names(coldata)) {
		stop("****** 'group' cannot be in the names of the description file, please change it ******")
	}
	all(rownames(coldata) == colnames(counts))
	return(coldata)
}

#READ gene desc file
makeDesc <- function(desc_gene) {
	desc<-read.table(file=desc_gene, sep="\t", header=F)
	names(desc)<-c("gene.id", "gene.name", "gene.type")
	desc <- desc[-grep("gene_id", desc$gene.id), ]
	return(desc)
}

filterDDS <- function(dds, mincount) {
	dds <- dds[ rowSums(counts(dds)) > mincount, ]
	dds <- DESeq(dds)
	return(dds)
}

# functions
ddsResult <- function(dds, desc, assay_field_raw, cond1, cond2 ) {
    assay_field <- str_remove(assay_field_raw, '~ ')
	filename<-paste(assay_field, cond1, cond2, sep="_")
	filename<-paste(filename, ".csv", sep="")
	res<-results(dds, contrast=c(assay_field,cond1,cond2))
	resOrdered <- res[order(res$padj),]
	resOrdered$ids<-row.names(resOrdered)
	resMerged<-merge(as.data.frame(resOrdered),desc, all = TRUE, by.x="ids", by.y="gene.id", sort=FALSE)
	resMerged$ID<-NULL
	write.csv(resMerged, file=filename, row.names = FALSE)          
    return (resMerged)
}

makeVolcano <- function(results, assay_field_raw, treat, ctrl, pcut, l2fc, h, w) {
    assayfield <- str_remove(assay_field_raw, '~ ')

	res_for_volc<-results[, c(8,3,7)]
	volcano_filename<-paste(assayfield, treat, ctrl, sep="_")
	volcano_file<-paste(volcano_filename, "_volc.pdf", sep="")

	subtitle<-paste(treat, ctrl, sep=" vs ")

	myplot<-EnhancedVolcano(res_for_volc,
    title = "Volcano plot",
    subtitle = subtitle,
    lab = res_for_volc$gene.name,
    x = 'log2FoldChange',
    pCutoff = pcut,
    FCcutoff = l2fc,
    y = 'padj')

	pdf(volcano_file, width = w, height = h)
    print(myplot)
	dev.off()
}


makeVST <- function(dds, blindval = FALSE) {
    vsd <- tryCatch( 
    {
      vst(dds, blind=blindval)
    },
      error = function(e) {
          varianceStabilizingTransformation(dds, blind=blindval)
      }
   )
   return(vsd)
}

sample_clustering <- function(vsd) {

	sampleDists  <- cor(assay(vsd), method = "spearman")
	sampleDists 

	sampleDistMatrix <- as.matrix(sampleDists )
	sampleDistMatrix

	return(sampleDistMatrix)
}


makeHeatmapDE <- function(results, dds, assay_field_raw, treat, ctrl, pcut, l2fc, h, w, maxnum=100000, sample_list=NULL, genes_file=NULL) {

    assayfield <- str_remove(assay_field_raw, '~ ')
    sig_res <- subset(results, padj <= pcut)
    if (maxnum < 100000) {
        sig_res<-sig_res[1:maxnum, ]
    }

    sig_ids <- sig_res$ids 

    # varianceStabilizingTransformation or vst?
    vsd<-assay(makeVST(dds, FALSE))
    sig_vsd_ids <- vsd[sig_ids, ]

    conv<-data.frame("ids"=results["ids"], "gene.name"=results["gene.name"])
    sig_vsd<-merge(as.data.frame(sig_vsd_ids),conv, all = FALSE, by.y="ids", by.x="row.names", sort=FALSE)
    sig_vsd["Row.names"]<-NULL
    row.names(sig_vsd)<-deduplicateIDs(sig_vsd$gene.name)
    sig_vsd["gene.name"]<-NULL

    if (!is.null(sample_list)) {
        samples<-str_split(sample_list, ",", simplify=TRUE)
        sig_vsd <- sig_vsd[ , samples]
    }

    if (!is.null(genes_file)) {
        mygenes<-trimws(readLines(genes_file))
	    mygenes <- mygenes[mygenes != ""]		
        sig_vsd <- sig_vsd[intersect(mygenes, rownames(sig_vsd)), ]
	    sig_size <- dim(sig_vsd)
	    print(paste0("we will analyze only ", sig_size[1], " genes"))
    } 

    heatmap_filename<-paste(assayfield, treat, ctrl, sep="_")
    heatmap_file<-paste(heatmap_filename, "_heat.pdf", sep="")

    subtitle<-paste("Heatmap of", treat, "vs",  ctrl, "( padj <=", pcut, ")",  sep=" ")
    pheatmap(sig_vsd, main=subtitle, cluster_rows=TRUE, filename = heatmap_file, show_rownames=TRUE, cluster_cols=TRUE, height=h, width=w)
}



deduplicateIDs <- function(vec) {
    dedup = ave(vec, vec, FUN =  \(x) if(length(x) > 1) paste(x, seq_along(x), sep = "__") else x)
    return(dedup)
}


plotPCA_plus <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, pcnum=6) {
    rv <- rowVars(assay(object), useNames=TRUE)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup,  drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else { colData(object)[[intgroup]] }
    d <- data.frame(group = group, intgroup.df, name = colnames(object))
	for (i in 1:pcnum) {
	d[[paste0("PC", i)]] <- pca$x[, i]
	}
	attr(d, "percentVar") <- percentVar[1:pcnum]
	pv <- round(100*attr(d, "percentVar"))
	PC_names <- c(paste("PC", 1:pcnum, sep=""))
	pv_df<-data.frame(PC_names, pv)

    if (returnData) {
 	    rots<-pca$rotation
 	    mylist<-list("data"=d,"rotation"=rots,"variance"=pv_df)
    	return(mylist)
	}
	

}

makePCA <- function(vsd, pca_fields, labels="", all_fields, PCA=1, PCB=2) {
    
    filename<-paste("PCA", pca_fields, sep="_")    
    filename<-paste(filename, "pdf", sep=".")    

	pdf(filename)
    PCA_lab<-paste("PC", PCA, sep="")
    PCB_lab<-paste("PC", PCB, sep="")
	cond=sym(pca_fields)
	PCA_labs=sym(PCA_lab)
	PCB_labs=sym(PCB_lab)
    pcaData_plus <- plotPCA_plus(vsd, intgroup = all_fields, returnData=TRUE)
    pcaData <- pcaData_plus[["data"]]
    rot.df<-as.data.frame(pcaData_plus[["rotation"]])
    rot.sel<-rot.df[, c(PCA_lab, PCB_lab)]
	percentVar <- round(100 * attr(pcaData, "percentVar"))
	myplot <- ggplot(pcaData, aes(!!PCA_labs, !!PCB_labs, color=!!cond))
	if (labels!="") {
		lab=sym(labels)
	    myplot <- ggplot(pcaData, aes(!!PCA_labs, !!PCB_labs, color=!!cond, label=!!lab)) +
	    geom_text(alpha = 0.8, size = 4, show.legend = FALSE)
	}
	myplot <- myplot + geom_point(size=2, alpha = 0.5) +
	xlab(paste0(PCA_lab, ": ",percentVar[PCA],"% variance")) +
	ylab(paste0(PCB_lab, ": ",percentVar[PCB],"% variance")) +
	scale_shape_identity() +
	scale_x_reverse() + theme_classic()
	print(myplot)
	dev.off()
	return(rot.sel)
}

printCounts <- function(dds, desc, ofile=NULL) {
    if (!is.null(ofile)) {
       vst_genes_f<-paste(ofile, "vst.genes", sep="_")
       norm_genes_f<-paste(ofile, "norm_counts.genes", sep="_")
       raw_genes_f<-paste(ofile, "raw_counts.genes", sep="_")
    } else {
       vst_genes_f<-"vst.genes"
       norm_genes_f<-"norm_counts.genes"
       raw_genes_f<-"raw_counts.genes"
    }

    vsd = makeVST(dds) 

    resMerged<-merge(as.data.frame(assay(vsd)), desc, all.x=F, by.x=0, by.y="gene.id", sort=FALSE)
    resMerged$ID<-NULL
    write.csv(resMerged, file=vst_genes_f, row.names = FALSE)          
    
    norcounts<-counts(dds, normalized=TRUE)
    rawcounts<-counts(dds, normalized=FALSE)
    resMerged2<-merge(as.data.frame(norcounts) ,desc, all.x=F, by.x=0, by.y="gene.id", sort=FALSE)

    resMerged2$ID<-NULL
    write.csv(resMerged2, file=norm_genes_f, row.names = FALSE)

    resMerged3<-merge(as.data.frame(rawcounts) ,desc, all.x=F, by.x=0, by.y="gene.id", sort=FALSE)
    resMerged3$ID<-NULL
 
    write.csv(resMerged3, file=raw_genes_f, row.names = FALSE)

    return(resMerged2)

}

printContrGenes <- function(rot_sel, desc) {

	resMerged<-merge(rot_sel, desc, all.x=F, by.x=0, by.y="gene.id", sort=FALSE)
	resOrdered <- resMerged[order(resMerged[, 2], decreasing = TRUE), ]
	write.csv(resOrdered, file="pca.genes", row.names = FALSE)          

}

## Function to uppercase gene names and remove any symbol except letters and numbers. 


normalize_gene_name <- function(x) {
  gsub("[^A-Z0-9]", "", toupper(x))
}


## Given a gene list in pca params, looks for the count in norm_counts, and write a yaml file with the counts grouped by gene and condition. This info will be given to multiqc to create gene boxplot for each group.

gene_expression_yaml <- function(norm_counts, genes, desc, intgroup = "condition") {
  
  ## Substitute the "." for "-" to match norm count colnames with desc$file names
  
  cols <- 2:(ncol(norm_counts) - 2)
  colnames(norm_counts)[cols] <- gsub("\\.", "-", colnames(norm_counts)[cols])
  

  
  ## ---- Identify columns ----
  norm_counts <- norm_counts[,-length(norm_counts)] ##Removing gene type column 
  gene_name_col <- colnames(norm_counts)[length(norm_counts)] 
  gene_id_col <- colnames(norm_counts)[1]
  
  ## Normalization of gene name column to match the gene list
  
  norm_counts$gene_name_norm <- normalize_gene_name(
    norm_counts[[gene_name_col]]
  )
  
  ## ---- Sanity checks ----
  if(!(intgroup %in% colnames(desc)))
    stop(paste("intgroup Column:",intgroup,"not found in your desc.txt"))
  
  ## Mapping interested genes in norm_counts
  
  #expr_filt <- norm_counts[norm_counts[[gene_name_col]] %in% genes, ]
  expr_filt <- norm_counts[
    norm_counts$gene_name_norm %in% genes,
  ]
  
  if (nrow(expr_filt) == 0)
    stop("genes not found in norm_counts")
  
  ## ---- Identify sample columns ----
  sample_cols <- setdiff(
    colnames(expr_filt),
    c(gene_id_col, gene_name_col,"gene_name_norm")
  )
  
  ## ---- Build file → condition lookup ----
  condition_map <- setNames(desc[[intgroup]], desc$file)
  
  ## ---- Build YAML structure ----
  yaml_list <- list()
  
  for (i in seq_len(nrow(expr_filt))) {
    
    gene <- expr_filt[[gene_name_col]][i]
    gene_entry <- list()
    
    for (col in sample_cols) {
      val <- expr_filt[i, col]
      
      if (!is.na(val) && val != 0) {
        condition <- condition_map[[col]]
        gene_entry[[condition]] <- c(gene_entry[[condition]], val)
      }
    }
    
    yaml_list[[gene]] <- gene_entry
  }
  
  ## ---- Write YAML ----
  yaml::write_yaml(yaml_list, "normalized_gene_counts_select.yaml")
  
  invisible(yaml_list)

  
}

genes_boxplot <- function(norm_counts,groups, genes_desc, genes, title) {



  gene_desc$gene_name_norm <- normalize_gene_name(gene_desc[,2])

  genes_match <- gene_desc[gene_desc$gene_name_norm %in% genes, ]

  if (nrow(genes_match) < 1) {
    print(paste0("Gene names: ",genes,"NOT found in gene_desc.txt"))}
  else {
  print(paste0("Gene name: ",genes_match$gene_name," found in gene_desc.txt"))

  #groups <- colData(vsd)[[1]]

  ## selecting vst counts for specific genes
  a <- vst_mat[rownames(vst_mat) %in% genes_match[[1]],]

  ab <- merge(a, genes_match, by.x=0, by.y="gene_id")
  rownames(ab) <- paste0(ab$gene_name,"_(",ab$Row.names,")")

  ## removing gene type, and gene name
  ab <- subset(ab, select=-c(Row.names,gene_name,gene_type,gene_name_norm))

  gene_cols <- rownames(ab)

  d <- t(ab)
  b <- colData(vsd)
  df <- merge(d,b, by = 0)

  ## Adding a new column with a short name of the sample for point labels
  df$sample_short <- sapply(
    strsplit(df$Row.names, "_"),
    function(x) paste(head(x, 2), collapse = "_")
  )

  ## Transforming table to have all the information all expression values
  df_long <- df |>
    pivot_longer(
      cols = all_of(gene_cols),
      names_to = "gene",
      values_to = "expression"
    )

  ## Creating a plot for each gene, wir one boxplot for condition
  ggplot(df_long, aes(x = condition, y = expression, fill = condition)) +
    geom_boxplot(width = 0.3 ,alpha = 0.20) +
    scale_fill_manual(values = c("black","#1C9BCD", "darkred","darkgreen", "#E69F00", "blue", "magenta")) +
    geom_jitter( shape = 21,              # supports fill + border
                 color = "black",         # border
                 alpha = 0.9,             # fill transparency
                 size = 2,
                 width = 0.2) +
    geom_text_repel(
      aes(label = sample_short),
      size = 2.5,
      show.legend = FALSE
    ) +
    facet_wrap(~ gene, nrow = 2, scales = "free_y") +
    labs(
      x = "Group",
      y = title,
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
  ggsave(paste(title,"genes_boxplot.png",sep="_"), width = 600,
         height =800 ,
         units = "px",
         dpi = 76)
  }
}

