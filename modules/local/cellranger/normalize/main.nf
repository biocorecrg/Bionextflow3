process CELLRANGER_NORM {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/nf-core/cellranger:9.0.1"

    input:
    tuple val(meta), path(pairs, name: "my-dir/*")

    output:
    tuple val(meta), path("*.gz"), emit: fastq


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """    
        ln -s ${pairs[0]} ${prefix}_S1_L001_R1_001.fastq.gz
        ln -s ${pairs[1]} ${prefix}_S1_L001_R2_001.fastq.gz 
    """
}
