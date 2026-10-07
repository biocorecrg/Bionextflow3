process CELLRANGER_COUNT {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/cellranger:10.0.0"

    input:
    tuple val(meta), path(pairs)
    path  index
    path  feature_reference

    output:
    tuple val(meta), path("**/outs", type: 'dir')   , emit: outs
    tuple val("${task.process}"), val('cellranger'), eval('cellranger --version | cut -d "-" -f 2'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def feature_ref_arg = feature_reference ? "--feature-ref=${feature_reference}" : ""
    """
    
    if [ ! -f ${prefix}_S1_L001_R1_001.fastq.gz ] && [ ! -f ${prefix}_S1_L001_R2_001.fastq.gz ]; then
        ln -s ${pairs[0]} ${prefix}_S1_L001_R1_001.fastq.gz
        ln -s ${pairs[1]} ${prefix}_S1_L001_R2_001.fastq.gz
    fi
        
    cellranger count ${args} --id=${prefix} \
                   --transcriptome=${index} \
                   --fastqs=./ \
                   --sample=${prefix} \
                   ${feature_ref_arg} \
                   --localcores=${task.cpus} \
                   --localmem=${task.memory.toGiga()} 
                   
    
   """
}
