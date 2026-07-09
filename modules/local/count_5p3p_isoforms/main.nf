process COUNT_5P3P_ISOFORMS {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py310h4a09ff2_1'
        : 'biocontainers/pysam:0.24.0--py310h4a09ff2_1'}"

    input:
    tuple val(meta), path(bams)
    tuple val(meta2), path(known_mirnas_gff3)

    script:
    def args = task.ext.args ?: ''
    """
    count_5p3p_isoforms.py \
        --gff3 ${known_mirnas_gff3} \
        --bams ${bams} \
        --output isoform_prop.tsv \
        --output-counts isoform_reads.tsv \
        ${args}
    """

    output:
    tuple val(meta), path("isoform_prop.tsv"), path("isoform_reads.tsv"), emit: isoform_counts
}
