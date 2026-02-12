process ISOQUANT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.10.0--hdfd78af_0':
        'biocontainers/isoquant:3.10.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(index)
    path(genome)
    path(annotation)

  
    output:
    tuple val(meta), path("OUT/*"), emit: out
    tuple val("${task.process}"), val("isoquant"), eval("isoquant.py --version | cut -d ' ' -f 2"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genedb_cmd = annotation ? "--genedb ${annotation}" : ''

    """
    isoquant.py --threads ${task.cpus} \\
    --reference ${genome} \\
    ${genedb_cmd} \\
    --bam ${bam} \\
    ${args} \\
    -o ./


    """


}
