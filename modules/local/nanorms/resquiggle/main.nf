process NANORMS_RESQUIGGLE {
    tag "$meta.id"
    label 'process_middle'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/lpryszcz/nanorms4:v4.0a2' :
        'docker.io/lpryszcz/nanorms4:v4.0a2' }"

    input:
    tuple val(meta), path(alnfile), path(alnindex)
    tuple val(meta2), path(pods)

    output:
    tuple val(meta), path("resquiggle/*.bam"),            emit: res_bam
    tuple val(meta), path("resquiggle/*.csi"),            emit: index
    tuple val("${task.process}"), val("bam2featureKS"), eval("/opt/app/src/resquiggle.py --version"), topic: versions, emit: versions_bam2featureks

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-vs-${meta2.id}"

    """
        /opt/app/src/resquiggle.py ${args} -t ${task.cpus} -k /opt/app/kmer_models/rna004.txt -o ./resquiggle -i ./ -b ${alnfile}
    """


}
