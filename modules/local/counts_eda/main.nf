process COUNTS_EDA {
    tag "$meta.id"
    container 'ghcr.io/biocorecrg/biocorecrg_deseq:latest'

    input:
    tuple val(meta), path(input_files)
    path(desc)
    path(dgenes)
    val(wtype)
    val(extrapars)
	
    output:
    // Main analysis outputs
    path("*_data.tsv"), emit: data
    path("*_variance.tsv"), emit: variance
    path("norm_counts.genes"), emit: norm_counts
    path("raw_counts.genes"), emit: raw_counts
    path("*Sample_Clustering_matrix.tsv"), emit: sample_clustering
    
    // Plot outputs (optional)
    path("_vst*boxplot.png"), emit: vst_genes_boxplot, optional: true
    path("_log2_deseq_norm*boxplot.png"), emit: norm_genes_boxplot, optional: true
    path("batch_vst*boxplot.png"), emit: batch_vst_genes_boxplot, optional: true
    path("batch_log2_deseq_norm*boxplot.png"), emit: batch_norm_genes_boxplot, optional: true
    
    // MultiQC-specific outputs (separate channel)
    path("*.{png,,log,json}"), emit: multiqc_files, optional: true
    
    // Versions (corrected syntax)
    path("versions.yml"), emit: versions, topic: versions


    script:
    """
    multiqc_pca.R -type ${wtype} -dgenes ${dgenes} -desc ${desc} ${extrapars}
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
        r: \$(R --version | head -n1 | cut -d' ' -f3)
    END_VERSIONS
    """
}