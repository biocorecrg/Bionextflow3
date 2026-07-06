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
    def software = meta.software ?: ''
    def filename = task.ext.args ?: 'supported_kmers'
    def out_name = software ? "${prefix}_${software}_${filename}" : "${prefix}_${filename}"

    """
	concatenate_bedrmod.R ${bedrmods_files} ${out_name}
    """

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}${meta.software ? '_' + meta.software : ''}_${task.ext.args ?: 'supported_kmers'}.*"), emit: results, optional: true
    tuple val("${task.process}"), val("merge_bedrmods"), eval("echo 1.0"), topic: versions
}
