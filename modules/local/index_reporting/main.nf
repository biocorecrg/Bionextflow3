process INDEX_REPORTING {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.11'
        : 'quay.io/biocontainers/python:3.11'}"

    input:
    tuple val(meta), path(llm_response_json), path(input_json)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_cmd = (input_json && input_json.name != 'NO_FILE') ? "--input-json ${input_json}" : ""
    """
    index_reporting.py \\
        --response-json ${llm_response_json} \\
        --output-mqc ${prefix}_index_reporting_mqc.html \\
        --sample-name "${prefix}" \\
        ${input_cmd} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "<!-- id: 'llm_evaluation' section_name: 'LLM evaluation' --><div>${prefix}: PASS</div>" > ${prefix}_index_reporting_mqc.html
    """

    output:
    tuple val(meta), path("${prefix}_index_reporting_mqc.html"), emit: mqc_html
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), topic: versions, emit: versions
}
