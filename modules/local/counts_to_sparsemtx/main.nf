process COUNTS_TO_SPARSEMTX {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/biocorecrg/bionextflow3/counts_to_sparsemtx:edge'

    input:
    tuple val(meta), path(counts, stageAs: "counts/*")
    path desc

    script:
    def args = task.ext.args ?: ''

    """
    counts_to_sparsemtx.py ${args} -i ./counts -o ${meta.id}_matrix/ -d ${desc}
    """

    output:
    tuple val(meta), path("*_matrix"), emit: matrix
}
