process LLM_SERVICE {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(json_file)
    path system_prompt
    path user_prompt

    output:
    tuple val(meta), path("${prefix}_llm_report.md")      , optional: true, emit: report
    tuple val(meta), path("${prefix}_llm_response.json")  , optional: true, emit: json
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), optional: true, topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def system_cmd = system_prompt ? "--system-prompt ${system_prompt}" : ""
    def user_cmd = user_prompt ? "--user-prompt ${user_prompt}" : ""
    """
    echo "Querying LLM web service..."
    query_llm_service.py \\
        --json-file ${json_file} \\
        --output-md ${prefix}_llm_report.md \\
        --output-json ${prefix}_llm_response.json \\
        ${system_cmd} \\
        ${user_cmd} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "# LLM Stub Report for ${prefix}" > ${prefix}_llm_report.md
    echo '{"choices":[{"message":{"content":"Stub report"}}]}' > ${prefix}_llm_response.json
    """
}
