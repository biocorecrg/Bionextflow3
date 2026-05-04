process SPACEMAKE_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/biocorecrg/spacemake:0.9.1' :
        'docker.io/biocorecrg/spacemake:0.9.1' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(comb_csv)
    tuple val(meta3), path(puck)
    tuple val(meta4), path(config)
    tuple val(meta5), path(barcode_folder)
    tuple val(meta6), path(genome)
    tuple val(meta7), path(annotation)

    output:
    tuple val("${task.process}"), val('spacemake'), eval('spacemake --version --version'), topic: versions, emit: versions
    tuple val(meta), path("project_df.csv") , emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:

    """
	spacemake run --cores ${task.cpus}
	
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fc_tiles

    """
}
