/*
*  MTX merge 
*/

process MTX_MERGE {
    tag "$meta.id"
    label 'process_high_memory'

    container 'docker://biocorecrg/sc_benchmark:0.2'

    input:
    tuple val(meta), path("*")

    output:
    tuple val(meta), path("*.h5seurat"), emit: merged_h5,    optional:true
    tuple val(meta), path("*.h5ad"),     emit: merged_h5ad,  optional:true

    script:
    def args           = task.ext.args ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"

    """
    merge_mtx.R ${args} -d ./ -o ${meta.id}
    """
} 	

