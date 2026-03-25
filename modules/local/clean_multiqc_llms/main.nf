process CLEAN_LLMS_PROMPT {
    tag "${meta.id}"
    container "biocontainers/python:3.10"

    label 'process_single'

    input:
    tuple val(meta), path(llms_file)
    val add_text

    script:

    """
    head -n 14 ${llms_file} > llms.system.txt
    tail -n +35 ${llms_file} | grep -v "mqc-custom-content-image" >> llms-full-cleaned.txt
    echo "${add_text}" >> llms.system.txt
    """

    output:
    tuple val(meta), path("llms-full-cleaned.txt"), emit: clean_llms
    tuple val(meta), path("llms.system.txt"), optional: true, emit: system_propmt
}
