process LLAMA_RUN_PYTHON {
    tag "${meta.id}"
    label 'process_gpu'

    //container "docker.io/biocorecrg/llama:0.1"
    container "quay.io/nf-core/llama-cpp-python:0.1.9"

    input:
    tuple val(meta), path(system), path(prompt_file)
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
        --output output.md \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch output.txt
    """

    output:
    tuple val(meta), path("output.md"), emit: output
}
