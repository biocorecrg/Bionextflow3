process DORADO_BASECALL {
    tag "$meta.id"
    label 'gpu'

    container "ontresearch/dorado:sha00aa724a69ddc5f47d82bd413039f912fdaf4e77"

    input:
    tuple val(meta), path(pod5)
    path  models
    
    output:
    tuple val(idfile), path("*.bam")        , emit: basecalled_bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    args = task.ext.args ?: ''
 
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    dorado basecaller  --models-directory ./ ${args} ./ > ${prefix}.bam                 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -n1 | cut -d "-" -f 2) 
END_VERSIONS
   """
}
