process DORADO_DOWNLOAD_MODEL {
    tag "$meta.id"
    label 'process_high'

    container "ontresearch/dorado:sha00aa724a69ddc5f47d82bd413039f912fdaf4e77"

    input:
    tuple val(meta), path(pod5)
    path models
    
    output:
    path models                 , type:'dir', emit: modelfolder
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    args = task.ext.args ?: ''
 
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    dorado basecaller ${args} --max-reads 1 --models-directory \$PWD/${models} ./ > test.bam; 
     
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -n1 | cut -d "-" -f 2) 
END_VERSIONS
   """
}

