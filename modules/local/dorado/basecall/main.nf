process DORADO_BASECALL {
    tag "$meta.id"
    label 'gpu'

    container "ontresearch/dorado:sha00aa724a69ddc5f47d82bd413039f912fdaf4e77"

    input:
    tuple val(meta), path(pod5)
    path  models
    
    output:
    tuple val(meta), path("*.bam")          , emit: basecalled_bam
    tuple val("${task.process}"), val("dorado"), eval("dorado --version 2>&1 | head -n1 | cut -d '-' -f 2 | cut -d '+' -f 1"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    def args = task.ext.args ?: ''
 
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    dorado basecaller  --models-directory ./${models} ${args} ./ > ${prefix}.bam                 
    

    """
}
