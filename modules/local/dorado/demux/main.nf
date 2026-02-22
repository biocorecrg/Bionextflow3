process DORADO_DEMUX {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"

    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.bam")          , emit: basecalled_bam
    tuple val("${task.process}"), val("dorado"), eval("dorado --version 2>&1 | head -n1 | cut -d '-' -f 2 | cut -d '+' -f 1"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    def args = task.ext.args ?: ''
 
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dorado demux --emit-summary $args --threads ${task.cpus} --output-dir ./ ${bam}
    gzip *_summary.txt
    # Avoid the folder structure from dorado
    mv */*/*/bam_pass/*/*.bam .


    """
}
