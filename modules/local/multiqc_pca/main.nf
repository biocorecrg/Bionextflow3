process MULTIQC_PCA {
    
    container 'fabianandrade/deseq2:2026_1'

    input:
    tuple val(meta), path(input_files)
    path(desc)
    path(dgenes)
    val(wtype)
    val(extrapars)
	
    output:
    path("*_data.tsv"), emit: data
    path("*_variance.tsv"), emit: variance
    path("norm_counts.genes"), emit: norm_counts
    path("raw_counts.genes"), emit: raw_counts
    path("Sample_Clustering_matrix.tsv"), emit: sample_clustering
    path("vst*boxplot.png"), emit: vst_genes_boxplot, optional: true
    path("log2_deseq_norm*boxplot.png"), emit: norm_genes_boxplot, optional:true

	
	
    script:
    """
    multiqc_pca.R -type ${wtype} -dgenes ${dgenes} -desc ${desc} ${extrapars} 
    """
}

