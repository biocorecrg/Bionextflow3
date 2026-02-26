process EPINANO_VARIANTS {
    tag "$meta.id"
    label 'process_middle'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/biocorecrg/mopepinano:1.2.5' :
        'docker.io/biocorecrg/mopepinano:1.2.5' }"

    input:
    tuple val(meta), path(alnfile), path(alnindex)
    path(reference)
    path(faiidx)
    path(dict_index)

    output:
    tuple val(meta), path("*.per.site.csv.gz"), emit: per_site_vars,       optional: true
    tuple val("${task.process}"), val('epinano'), val('1.2.5'), topic: versions, emit: versions_epinano
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export XDG_CACHE_HOME="./tmp"
	Epinano_Variants_CRAM.py -c ${task.cpus} -r ${reference} -b ${alnfile}
        shopt -s nullglob
	for i in *.csv; do gzip \$i; done
    """


}
