process CELLRANGER_COUNT {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/nf-core/cellranger:9.0.1"

    input:
    tuple val(meta), path(pairs)
    path  index

    output:
    tuple val(meta), path("**/outs/**"), emit: outs
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    args = task.ext.args ?: ''
 
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    if [ ! -f ${prefix}_S1_L001_R1_001.fastq.gz ] && [ ! -f ${prefix}_S1_L001_R2_001.fastq.gz ]; then
        ln -s ${pairs[0]} ${prefix}_S1_L001_R1_001.fastq.gz
        ln -s ${pairs[1]} ${prefix}_S1_L001_R2_001.fastq.gz
    fi
        
	cellranger count ${args} --id=${prefix} \
                   --transcriptome=${index} \
                   --fastqs=./ \
                   --sample=${prefix} \
                   --localcores=${task.cpus} \
                   --localmem=${task.memory.toGiga()} 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
   """
}
