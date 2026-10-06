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
    tuple val(meta), path("${prefix}_demux_llm_mqc.html") , optional: true, emit: mqc_html
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), optional: true, topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // v2: formatted MultiQC table output
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def system_cmd = system_prompt ? "--system-prompt ${system_prompt}" : ""
    def user_cmd = user_prompt ? "--user-prompt ${user_prompt}" : ""
    """
    echo "Diagnosing demultiplexing stats with LLM..."
    query_llm_service.py \\
        --json-file ${json_file} \\
        --output-md ${prefix}_llm_report.md \\
        --output-json ${prefix}_llm_response.json \\
        --output-mqc ${prefix}_demux_llm_mqc.html \\
        --sample-name "${prefix}" \\
        ${system_cmd} \\
        ${user_cmd} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "# LLM Stub Report for ${prefix}" > ${prefix}_llm_report.md
    echo '{"choices":[{"message":{"content":"Stub report"}}]}' > ${prefix}_llm_response.json
    echo "<!-- id: 'demultiplexing_llm_evaluation' section_name: 'Demultiplexing LLM evaluation' --><div>${prefix}: PASS</div>" > ${prefix}_demux_llm_mqc.html
    """
}
