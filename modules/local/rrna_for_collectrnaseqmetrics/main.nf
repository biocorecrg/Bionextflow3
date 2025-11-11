// custom
process RRNA_FOR_COLLECTRNASEQMETRICS {
    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocorecrg/almalinux-perlbrew-pyenv3':
        'biocorecrg/almalinux-perlbrew-pyenv3' }"

    input:
    tuple val(meta), path (faidx)
    tuple val(meta2), path (gtf)

    output:
    tuple val(meta), path("rRNA.list")  , emit: rRRNAs
    
    script:
    """
    rrna_for_collect.sh ${faidx} ${gtf}
    """
}
