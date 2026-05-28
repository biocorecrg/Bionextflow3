process DE_PCA {
    tag "${meta.id}"
    container 'ghcr.io/biocorecrg/biocorecrg_deseq:latest'

    input:
    tuple val(meta), path(input_files)
    path desc
    path dgenes
    val wtype
    val extrapars

    script:
    """
    counts_eda.R -type ${wtype} -dgenes ${dgenes} -desc ${desc} ${extrapars}
     """

    output:
    // Main analysis outputs
    path ("*_data.tsv"), emit: data
    path ("*_variance.tsv"), emit: variance
    path ("norm_counts.genes"), emit: norm_counts
    path ("raw_counts.genes"), emit: raw_counts
    path ("*Sample_Clustering_matrix.tsv"), emit: sample_clustering

    // Plot outputs (optional)
    path ("_vst*boxplot.png"), emit: vst_genes_boxplot, optional: true
    path ("_log2_deseq_norm*boxplot.png"), emit: norm_genes_boxplot, optional: true
    path ("batch_vst*boxplot.png"), emit: batch_vst_genes_boxplot, optional: true
    path ("batch_log2_deseq_norm*boxplot.png"), emit: batch_norm_genes_boxplot, optional: true

    // GSEA outputs (optional)
    path ("gsea*.{txt,cls}"), emit: gsea_files, optional: true

    // MultiQC-specific outputs (separate channel)
    path ("*.{png,tsv,genes,json}"), emit: multiqc_files, optional: true

    //Versions 
    tuple val("${task.process}"), val("deseq2"), eval("Rscript -e 'cat(as.character(packageVersion(\"DESeq2\")))'"), topic: versions
    tuple val("${task.process}"), val("ggplot2"), eval("Rscript -e 'cat(as.character(packageVersion(\"ggplot2\")))'"), topic: versions
    tuple val("${task.process}"), val("r"), eval("R --version | head -1 | sed 's/R version //;s/ (.*//'"), topic: versions

}
