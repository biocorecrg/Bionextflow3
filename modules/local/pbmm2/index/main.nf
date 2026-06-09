process PBMM2_INDEX {
    tag "${fasta}"
    label 'process_low'

    // Note: the versions here need to match the versions used in the container below
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pbmm2_samtools:188cdf4b42a3fc20'
        : 'community.wave.seqera.io/library/pbmm2_samtools:afa56306dfbe0957'}"

    input:
    tuple val(meta), path(fasta)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    pbmm2 \\
        index \\
        ${args} \\
        ${fasta} \\
        ${prefix}.mmi
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.mmi
    """

    output:
    tuple val(meta), path("*.mmi"), emit: index
    tuple val("${task.process}"), val("pbmm2"), eval("pbmm2 --version"), topic: versions, emit: versions
    tuple val("${task.process}"), val("samtools"), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
}
