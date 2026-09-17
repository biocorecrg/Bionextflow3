process DORADO_BASECALL {
    tag "$meta.id"
    label 'gpu'

    container "docker://nanoporetech/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f"

    input:
    tuple val(meta), path(pod5)
    path  models
    
    output:
    tuple val(meta), path("*.bam")          , optional:true , emit: basecalled_bam
    tuple val(meta), path("*.cram")         , optional:true , emit: basecalled_cram
    tuple val("${task.process}"), val("dorado"), eval("dorado --version 2>&1 | head -n1 | cut -d '-' -f 2 | cut -d '+' -f 1"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    def args = task.ext.args ?: ''
 
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outext = args.contains("--emit-cram") ? "cram" : "bam"

    """
    dorado basecaller  --models-directory ./${models} ${args} ./ > ${prefix}.${outext}              
    """
}
