process CLAIRS {
    tag "${meta.id}"
    label 'gpu'
    label 'process_high'

    container "docker.io/hkubal/clairs:latest"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(reference_index)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clairs \\
        --tumor_bam_fn=${tumor_bam} \\
        --normal_bam_fn=${normal_bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --platform=${platform} \\
        --output_dir=. \\
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
    tuple val("${task.process}"), val('clairs'), eval('run_clairs --version | sed "s/^ClairS v//"'), emit: versions, topic: versions
}
