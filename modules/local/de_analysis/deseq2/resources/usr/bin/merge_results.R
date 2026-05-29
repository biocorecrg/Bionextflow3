#!/usr/local/bin/Rscript --vanilla

source("/functions.R")

# Load common pars
parser <- getCommonPars()

#Define desired outputs:
#SPECIFIC FEATURES:
parser$add_argument("-input", "--Input_f", type="character", help="Path or matrix file to read counts")
parser$add_argument("-type", "--Input_type", type="character", help="Type of input: counts, salmon, star, matrix")
parser$add_argument("-contrast", "--Contrast", type="character", help="value for contast", default="")
parser$add_argument("-comp", "--comparisons", type="character", help="Comma separated list of comparisons, groups separated by dash: treat-control (e.g: -comp R2WT-WT40,R2KO-WT40,KO40-WT40")
parser$add_argument("-strand", "--Strand_id", type="integer", default=4, help="Strand to analyze, 2: unstranded, 3: forward, 4 reverse (ONLY FOR STAR) [default %(default)s]")
parser$add_argument("-dtrans", "--Desc_Tx", type="character", default="tx2gene.csv", help="Path to read tx2gene file (ONLY FOR SALMON) [default %(default)s]")



#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

if (args$Contrast == "") {
	contrast = str_remove(args$Assay_field, '~ ')
} else contrast = args$Contrast
 

if (args$Input_type == "counts") {
	out <- makeDDSFromCounts(args$Input_f, args$Desc_exp, args$Assay_field)
} else if (args$Input_type == "salmon") {
	out <- makeDDSFromSalmon(args$Input_f, args$Desc_exp, args$Assay_field, args$Desc_Tx)
} else if (args$Input_type == "star") {
	out <- makeDDSFromStar(args$Input_f, args$Desc_exp, args$Assay_field, args$Strand_id)
} else if (args$Input_type == "matrix") {
	out <- makeDDSFromMatrix(args$Input_f, args$Desc_exp, args$Assay_field)
} else {
	print("Please define one correct input type!")
    stop()
}

dds <- filterDDS(out$dds, args$min_count)

# Get results
desc <- makeDesc(args$Desc_genes)

topall <- data.frame(ids = desc$gene.id)

## metadata

metadata <- read.csv(args$Desc_exp, sep="\t", header=T)

# Get results for all comparisons combinations if not provided by the user

if (length(args$comparisons) == 0) {
    groups = factor(metadata[[contrast]])
    ## replacing "-" for "_" in group names to avoid issues with column names in the final table
    levels(groups) <- gsub("-", "_", levels(groups))
    levels(groups) <- unique(groups)
    comparisons <- c()
    for (i in 1:(length(levels(groups)) - 1)) {
        for (j in (i + 1):length(levels(groups))) {
            comparisons <- c(comparisons, paste0(levels(groups)[j], "-", levels(groups)[i]))
        }
    
    }
    print(paste0("Creating results for the comparison: ",comparisons))
} else {
    comparisons <- strsplit(args$comparisons, ",")[[1]]
    print(paste0("Creating results for the comparison: ",comparisons))
}


## Creating results table for each comparison and merging them into one final table

for (i in 1:length(comparisons)) {
    comparison <- comparisons[i]
    groups <- strsplit(comparison, "-")[[1]]
    treat_group <- groups[1]
    control_group <- groups[2]
    results <- ddsResult(dds, desc, contrast, treat_group, control_group)
    rownames(results) <- results$ids

    log2_comp <- paste0("log2FoldChange_", treat_group, "_vs_", control_group)
    pvalue_comp <- paste0("pvalue_", treat_group, "_vs_", control_group)
    padj_comp <- paste0("padj_", treat_group, "_vs_", control_group)
 
    results$gene.type <- NULL
    results$gene.id <- NULL
    results$gene.name <- NULL
    results$chromosome <- NULL
    results$start <- NULL
    results$end <- NULL
    results$strand <- NULL
    results$baseMean <- NULL
    results$lfcSE <- NULL
    results$stat <- NULL


    names(results)[names(results) == "log2FoldChange"] <- log2_comp
    names(results)[names(results) == "pvalue"] <- pvalue_comp
    names(results)[names(results) == "padj"] <- padj_comp

    topall <- merge(topall, results, by = "ids")
    
}

## Final table with annotattion and results for all comparisons

final_merge <- merge(desc, topall, by.x = "gene.id", by.y = "ids")
write.table(final_merge, file = "merged_results.csv", sep = ",", row.names = FALSE)