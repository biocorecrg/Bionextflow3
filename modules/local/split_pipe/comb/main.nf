process SPLITPIPE_COMB {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/biocorecrg/spipe:1.6.1'

    input:
    tuple val(meta),   path(sublib_folders)
    tuple val(meta2),  path(index)
    
    output:
    tuple val(meta), path("${meta.id}")                       , emit: out
    path  "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    split-pipe \
       --mode comb --sublibraries ${sublib_folders} \
       --output_dir ${meta.id}
  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split-pipe: \$(split-pipe --version | cut -d ' ' -f 2 | sed 's/v//g')
END_VERSIONS
    """


}


