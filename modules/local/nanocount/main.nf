process NANOCOUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocount:1.0.0.post6--pyhdfd78af_0':
        'biocontainers/nanocount:1.0.0.post6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(index)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    NanoCount \\
    -i ${input} \\
    -o ${prefix}.txt \\
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocount: \$(echo \$(NanoCount --version ))
    END_VERSIONS
    """


}
