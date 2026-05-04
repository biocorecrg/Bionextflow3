process FILTER_SHORTSTACK {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_filter.bam")                       , emit: filtered_bam
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed "s/samtools //"'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
    samtools view -@ ${task.cpus} -h ${bam} | awk '\$0 ~ /^@/ || \$0 ~ /XY:Z:U/' | samtools view -bS - > ${meta.id}_filter.bam

    """

}
