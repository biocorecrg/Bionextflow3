process CELLRANGER_MULTI {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/nf-core/cellranger:9.0.1"

    input:
    tuple val(meta), path(pairs), val(barcodes)
    path index

    output:
    tuple val(meta), path("**/outs/**"), emit: outs
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    args = task.ext.args ?: ''
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def list_ids = barcodes.join("\n")  

    """
    
    if [ ! -f ${prefix}_S1_L001_R1_001.fastq.gz ] && [ ! -f ${prefix}_S1_L001_R2_001.fastq.gz ]; then
        ln -s ${pairs[0]} ${prefix}_S1_L001_R1_001.fastq.gz
        ln -s ${pairs[1]} ${prefix}_S1_L001_R2_001.fastq.gz
    fi
    
    echo  "[gene-expression]" > multi_config.csv
    echo  "reference,\$PWD/${index}" >> multi_config.csv
    echo  "create-bam,false" >> multi_config.csv
    echo  "" >> multi_config.csv
    echo  "[libraries]" >> multi_config.csv
    echo  "fastq_id,fastqs,feature_types" >> multi_config.csv
    echo  "${prefix},\$PWD,Gene Expression" >> multi_config.csv
    echo  "" >> multi_config.csv
    echo  "[samples]" >> multi_config.csv
    echo  "sample_id,ocm_barcode_ids" >> multi_config.csv
    echo  "${list_ids}" >> multi_config.csv

    cellranger multi --id=${prefix}  ${args} --csv=multi_config.csv  \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()}     
                
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version | cut -d "-" -f 2) 
END_VERSIONS
   """
}
