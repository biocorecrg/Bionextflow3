process LLAMA_RUN_PYTHON {
    tag "${meta.id}"
    label 'process_gpu'

    //container "docker.io/biocorecrg/llama:0.1"
    container "quay.io/nf-core/llama-cpp-python:0.1.9"

    input:
    tuple val(meta), path(prompt_file), path(system)
    path model

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    llamacpp-python.py \
        --model ${model} \
        --system_prompt ${system} \
        --messages ${prompt_file} \
        --output ${prefix}.md \
        --html_output ${prefix}_mqc.html \
        --sample_name "${prefix}" \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.md
    touch ${prefix}_mqc.html
    """

    output:
    tuple val(meta), path("${prefix}.md"), emit: output
    tuple val(meta), path("${prefix}_mqc.html"), emit: html
}
