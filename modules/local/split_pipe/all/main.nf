process SPLITPIPE_ALL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'biocorecrg/spipe:1.3.1'

    input:
    tuple val(meta),   path(reads)
    tuple val(meta2),  path(index)
    
    output:
    tuple val(meta), path("${meta.id}")                       , emit: out
    path  "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def reads1 = []
    def reads2 = []
    meta.single_end ? [reads].flatten().each { reads1 << it } : reads.eachWithIndex { v, ix -> (ix & 1 ? reads2 : reads1) << v }
    def input_reads = meta.single_end ? "-r ${reads1.join(" ")}" : "--fq1 ${reads1.join(" ")} --fq2 ${reads2.join(" ")}"
  
    """
   split-pipe \
       --mode all ${args}  \
       --genome_dir ${index} \
       ${input_reads} \
       --output_dir ${prefix} \
       --nthreads ${task.cpus}

  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split-pipe: \$(split-pipe --version | cut -d ' ' -f 2 | sed 's/v//g')
    END_VERSIONS
    """

}
