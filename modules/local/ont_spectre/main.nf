process ONT_SPECTRE {
    tag { "${meta.id}" }
    label 'process_high'

    container "docker.io/ontresearch/spectre:sha42472d37a5a992c3ee27894a23dce5e2fff66d27"

    input:
    tuple val(meta), path("mosdepth/*"), path(reference), path(vcf)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    spectre CNVCaller \\
        ${args} \\
        --bin-size 1000 \\
        --coverage mosdepth \\
        --sample-id ${meta.id} \\
        --output-dir ./${meta.id} \\
        --reference ${reference} \\
        --snv ${vcf}
    """

    output:
    tuple val(meta), path("${meta.id}/"), emit: outdir
    tuple val(meta), path("${meta.id}/${meta.id}.vcf"), emit: vcf, optional: true
    tuple val("${task.process}"), val('spectre'), eval('pip show ont-spectre | grep Version | sed "s/Version: //"'), topic: versions, emit: versions
}
