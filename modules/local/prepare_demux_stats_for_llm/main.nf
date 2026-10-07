process PREPARE_DEMUX_STATS_FOR_LLM {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(stats_json)

    output:
    tuple val(meta), path("${prefix}_demux_stats_llm.json"), emit: json
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), optional: true, topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    prepare_demux_stats_for_llm.py \\
        --stats-json ${stats_json} \\
        --output-json ${prefix}_demux_stats_llm.json \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{}' > ${prefix}_demux_stats_llm.json
    """
}
