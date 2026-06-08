process ADAPTIVE_SAMPLINGQC {
    label 'process_low'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.10'
        : 'quay.io/biocontainers/python:3.10'}"

    input:
    tuple val(meta), path(as_decisions), path(sequencing_summary)

    script:
    args = task.ext.args ?: ''

    """
    AS_QC.py \
        ${args} \
        --as_decisions ${as_decisions} \
        --sequencing_summary ${sequencing_summary} \
        --output_name ${meta.id}
    """

    output:
    tuple val(meta), path("*-Accepted-Reads.txt"), emit: accepted_reads
    tuple val(meta), path("*_mqc.txt"), emit: stats
}
