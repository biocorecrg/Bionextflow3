process MERGE_BEDRMODS {
    tag "${meta.id}"
    label 'process_short_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker.io/biocorecrg/mop_consensus:0.1'
        : 'docker.io/biocorecrg/mop_consensus:0.1'}"

    input:
    tuple val(meta), path(bedrmods_files)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id

    """
	concatenate_bedrmod.R ${bedrmods_files} ${prefix}_supported_kmers.bedrmod
    """

    output:
    tuple val(meta), path("*_supported_kmers.bedrmod"), emit: results
    tuple val("${task.process}"), val("merge_bedrmods"), eval("echo 1.0"), topic: versions
}