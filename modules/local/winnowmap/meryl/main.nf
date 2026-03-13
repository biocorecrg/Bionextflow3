process WINNOWMAP_MERYL {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/samtools_winnowmap:5e97225b2810c151' :
        'community.wave.seqera.io/library/samtools_winnowmap:5e97225b2810c151' }"

    input:
    tuple val(meta), path(reference)
 
    output:
    tuple val(meta), path("repetitive_k.txt"),                      emit: repetitivek
    tuple val("${task.process}"), val("winnowmap"), eval("winnowmap --version"), topic: versions, emit: versions

 
    script:
    def args  = task.ext.args ?: ''

    """
    meryl count $args output merylDB ${reference}
  	meryl print greater-than distinct=0.9998 merylDB > repetitive_k.txt	
    """

    stub:
    """
    touch repetitive_k.txt	

    """
}
