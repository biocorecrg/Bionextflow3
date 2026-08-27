process CLAIR3_RNA {
    tag "${meta.id}"
    label 'gpu'
    label 'process_high'

    container "docker.io/hkubal/clair3-rna:latest"

    input:
    tuple val(meta), path(bam), path(bai), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(reference_index)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clair3_rna \\
        -B ${bam} \\
        -R ${reference} \\
        -o . \\
        -t ${task.cpus} \\
        -p ${platform} \\
        --output_prefix=${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val("${task.process}"), val('clair3-rna'), eval('run_clair3_rna --version | sed "s/^Clair3-RNA v//"'), emit: versions, topic: versions
}
