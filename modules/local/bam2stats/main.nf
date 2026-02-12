process BAM2STATS {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1a35167f7a491c7086c13835aaa74b39f1f43979:9254eac8981f615fb6c417fa44e77c3b44bc3abd-0' :
        'quay.io/biocontainers/mulled-v2-1a35167f7a491c7086c13835aaa74b39f1f43979:a7b00ff483a30f0a985d9e0d4da1f5762af68cd6-0' }"

    input:
    tuple val(meta), path(bamfile)

    output:
    tuple val(meta), path("*.stat"), emit: stats
    tuple val("${task.process}"), val("samtools"), eval("samtools --version | head -n1 | cut -d ' ' -f 2"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
        bam2stats.py ${bamfile} > ${prefix}.stat
    """

}
