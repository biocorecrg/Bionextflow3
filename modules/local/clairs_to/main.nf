process CLAIRS_TO {
    tag "${meta.id}"
    label 'gpu'
    label 'process_high'

    container "docker.io/hkubal/clairs-to:latest"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(reference_index)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clairs_to \\
        --tumor_bam_fn=${tumor_bam} \\
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
    echo "" | gzip > ${prefix}.snv.vcf.gz
    touch ${prefix}.snv.vcf.gz.tbi
    echo "" | gzip > ${prefix}.indel.vcf.gz
    touch ${prefix}.indel.vcf.gz.tbi
    """

    output:
    tuple val(meta), path("${prefix}.snv.vcf.gz"), emit: snv_vcf, optional: true
    tuple val(meta), path("${prefix}.snv.vcf.gz.tbi"), emit: snv_tbi, optional: true
    tuple val(meta), path("${prefix}.indel.vcf.gz"), emit: indel_vcf, optional: true
    tuple val(meta), path("${prefix}.indel.vcf.gz.tbi"), emit: indel_tbi, optional: true
    tuple val("${task.process}"), val('clairs-to'), eval('run_clairs_to --version | sed "s/^ClairS-TO v//"'), emit: versions, topic: versions
}
