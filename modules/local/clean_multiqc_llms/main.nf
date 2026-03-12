process CLEAN_LLMS_PROMPT {
    tag "${meta.id}"
    container "biocontainers/python:3.10"

    label 'process_single'


    input:
    tuple val(meta), path(llms_file)
    
    output:
    tuple val(meta), path("llms-full-cleaned.txt")
    
    script:
    """
    head -n 13 ${llms_file} > llms-full-cleaned.txt
    tail -n +35 ${llms_file} >> llms-full-cleaned.txt
    """
}
