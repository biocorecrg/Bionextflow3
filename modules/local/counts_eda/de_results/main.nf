process DE_RESULTS {
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
    merge_results.R -type ${wtype} -dgenes ${dgenes} -desc ${desc} ${extrapars}
    """

    output:
    // Individual results table for each comparisons.
    path ("*.csv"), emit: individual_results

    // Merged results table combining all comparisons. 
    path ("merged_results.csv"), emit: merged_results

}
