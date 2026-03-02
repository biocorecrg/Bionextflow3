process NANORMS_BAM2FEATUREKS {
    tag "$meta.id vs $meta2.id"
    label 'process_middle'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/lpryszcz/nanorms4:v4.0a2' :
        'docker.io/lpryszcz/nanorms4:v4.0a2' }"

    input:
    tuple val(meta), path(alnfileA), path(alnindexA)
    tuple val(meta2), path(alnfileB), path(alnindexB)

    output:
    tuple val(meta), path("*_baseQ.bed.gz"),            emit: bed
    tuple val("${task.process}"), val("bam2featureKS"), eval("/opt/app/src/bam2featureKS.py --version"), topic: versions, emit: versions_bam2featurehs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-vs-${meta2.id}"
    """
        /opt/app/src/bam2featureKS.py ${args}  -t ${task.cpus} -i ${alnfileA} ${alnfileB} -o ./${prefix}
    """


}
