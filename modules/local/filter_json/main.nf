process FILTER_JSON {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jq:1.7--h5bf99c6_0' :
        'quay.io/biocontainers/jq:1.7--h5bf99c6_0' }"

    input:
    tuple val(meta), path(json_file)
    val filter

    output:
    tuple val(meta), path("${prefix}_filtered.json"), emit: json
    tuple val("${task.process}"), val('jq'), eval('jq --version | sed "s/jq-//"'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    jq '${filter}' ${json_file} > ${prefix}_filtered.json
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{}' > ${prefix}_filtered.json
    """
}
