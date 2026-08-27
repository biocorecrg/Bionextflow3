process CLAIR_MOSAIC {
    tag "${meta.id}"
    label 'gpu'
    label 'process_high'

    container "docker.io/hkubal/clair-mosaic:latest"

    input:
    tuple val(meta), path(bam), path(bai), path(control_bam), path(control_bai), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(reference_index)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}."
    def control = control_bam ? "--control_bam_fn=${control_bam}" : ""
    """
    run_clair_mosaic \\
        --bam_fn=${bam} \\
        ${control} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --output_dir=. \\
        --platform=${platform} \\
        ${args}

    # Rename to add prefix
    for file in snv.vcf.gz snv.vcf.gz.tbi indel.vcf.gz indel.vcf.gz.tbi; do
        if [ -e "\$file" ]; then
            mv "\$file" "${prefix}\$file"
        fi
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}."
    """
    echo "" | gzip > ${prefix}snv.vcf.gz
    touch ${prefix}snv.vcf.gz.tbi
    echo "" | gzip > ${prefix}indel.vcf.gz
    touch ${prefix}indel.vcf.gz.tbi
    """

    output:
    tuple val(meta), path("${prefix}snv.vcf.gz"), emit: snv_vcf, optional: true
    tuple val(meta), path("${prefix}snv.vcf.gz.tbi"), emit: snv_tbi, optional: true
    tuple val(meta), path("${prefix}indel.vcf.gz"), emit: indel_vcf, optional: true
    tuple val(meta), path("${prefix}indel.vcf.gz.tbi"), emit: indel_tbi, optional: true
    tuple val("${task.process}"), val('clair-mosaic'), eval('run_clair_mosaic --version | sed "s/^Clair-Mosaic v//"'), emit: versions, topic: versions
}
