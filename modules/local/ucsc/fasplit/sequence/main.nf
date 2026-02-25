process UCSC_FASPLIT_SEQUENCE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-fasplit:482--h0b57e2e_0' :
        'biocontainers/ucsc-fasplit:482--h0b57e2e_0' }"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    faSplit sequence ${sequence} ${args} ${prefix}
    
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.num.fa
    """
}
