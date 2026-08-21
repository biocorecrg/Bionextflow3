process CELLRANGER_MULTI {
    tag "${meta.id}"
    label 'process_high'

    container "nf-core/cellranger:10.0.0"

    input:
    tuple val(meta), path(pairs), val(barcodes)
    path index
    path probes
    path feature_reference

    when:
    task.ext.when == null || task.ext.when

    script:

    args = (task.ext.args ?: 'create-bam,false').replaceAll(' ', '\n')

    def prefix = task.ext.prefix ?: "${meta.id}"
    def list_ids = barcodes ? barcodes.join("\n") : ""

    """
    
    if [ ! -f ${prefix}_S1_L001_R1_001.fastq.gz ] && [ ! -f ${prefix}_S1_L001_R2_001.fastq.gz ]; then
        ln -s ${pairs[0]} ${prefix}_S1_L001_R1_001.fastq.gz
        ln -s ${pairs[1]} ${prefix}_S1_L001_R2_001.fastq.gz
    fi
    
    echo  "[gene-expression]" > multi_config.csv
    if [ -n "${index}" ]; then
        echo  "reference,\$PWD/${index}" >> multi_config.csv
    fi
    if [ -n "${probes}" ]; then
        echo  "probe-set,\$PWD/${probes}" >> multi_config.csv
    fi
    echo  "${args}" >> multi_config.csv
    echo  "" >> multi_config.csv

    if [ -n "${feature_reference}" ]; then
        echo  "[feature]" >> multi_config.csv
        echo  "reference,\$PWD/${feature_reference}" >> multi_config.csv
        echo  "" >> multi_config.csv
    fi

    echo  "[libraries]" >> multi_config.csv
    echo  "fastq_id,fastqs,feature_types" >> multi_config.csv
    echo  "${prefix},\$PWD,Gene Expression" >> multi_config.csv
    echo  "" >> multi_config.csv

    if [ -n "${list_ids}" ]; then
        echo  "[samples]" >> multi_config.csv
        echo  "sample_id,ocm_barcode_ids" >> multi_config.csv
        echo  "${list_ids}" >> multi_config.csv
    fi

    cellranger multi --id=${prefix} --csv=multi_config.csv  \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()}     
                
   """

    output:
    tuple val(meta), path("**/outs/per_sample_outs/*", type: 'dir'), emit: outs
    tuple val("${task.process}"), val('cellranger'), eval('cellranger --version | cut -d "-" -f 2'), topic: versions, emit: versions
}
