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
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    isoquant.py --threads ${task.cpus} \\
    --reference ${genome} \\
    --genedb ${annotation} \\
    --bam ${bam} \\
    ${args}
    -o ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(echo \$(isoquant.py --version ))
    END_VERSIONS
    """


}
