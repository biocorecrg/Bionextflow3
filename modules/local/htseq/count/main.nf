process HTSEQ_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.3--py310ha14a713_0':
        'biocontainers/htseq:2.0.3--py310ha14a713_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.txt")                      , emit: txt
    tuple val(meta), path("*_anno.bam"),  optional: true, emit: bam
    tuple val("${task.process}"), val("htseq"), eval("htseq-count --version | sed 's/^.*htseq-count //; s/Using.*\$//'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_output = args.contains("-p bam") ? "-o ${prefix}_anno.bam" : ""

    """
    htseq-count \\
        ${input} \\
        ${gtf} \\
        ${args} \\
        ${bam_output} \\
        > ${prefix}.txt



    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt


    """
}
