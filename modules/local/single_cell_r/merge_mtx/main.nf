/*
*  MTX merge 
*/

process MTX_merge {
    tag "$meta.id"
    label 'process_high_memory'

    container 'biocorecrg/sc_benchmark:0.2'

    input:
    tuple val(meta), path("*")

    output:
    tuple val(meta), path("*.h5"), emit: h5

    script:
    args = task.ext.args ?: ''
 
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    merge_mtx.R -d ./ -o ${meta.id}
    """
} 	

