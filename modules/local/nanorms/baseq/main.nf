process NANORMS_BASEQ {
    tag "$meta.id vs $meta2.id"
    label 'process_middle'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/lpryszcz/nanorms4:v4.0a2' :
        'docker.io/lpryszcz/nanorms4:v4.0a2' }"

    input:
    tuple val(meta), path(alnfileA), path(alnindexA)
    tuple val(meta2), path(alnfileB), path(alnindexB)

    output:
    tuple val(meta), path("*_baseQ.bed.gz"),         emit: bed
    tuple val("${task.process}"), val('nanorms'), val('4.0.a2'), topic: versions, emit: versions_nanorms
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-vs-${meta2.id}"
    """
      	/opt/app/src/bam2baseQ.py ${args} -i ${alnfileA} ${alnfileB} -t4 | gzip -c > ${prefix}_baseQ.bed.gz
   """


}
