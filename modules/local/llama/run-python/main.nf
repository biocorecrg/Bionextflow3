process LLAMA_RUN_PYTHON {
    tag "$meta.id"
    label 'process_gpu'

    container "docker.io/biocorecrg/llama:0.1"
    //container "community.wave.seqera.io/library/llama-cpp-python_cuda:e081d1272ef7c140"
    //container "ghcr.io/abetlen/llama-cpp-python:latest"

    input:
    tuple val(meta), path(prompt_file)
    path(model)

    output:
    tuple val(meta), path("output.md"), emit: output
    //tuple val("${task.process}"), val("llama"), eval("llama: \$(python3 -c 'import llama_cpp; print(llama_cpp.__version__)')"), topic: versions, emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    llamacpp-python.py \
        --model ${model} \
        --messages ${prompt_file} \
        --output output.md \
        ${args} 
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch output.txt
    """
}
