process PBMM2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pbmm2_samtools:188cdf4b42a3fc20'
        : 'community.wave.seqera.io/library/pbmm2_samtools:afa56306dfbe0957'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(fasta)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_fmt = "bam"
    if (args2.contains('--output-fmt')) {
        def fmt_matcher = args2 =~ /--output-fmt\s+([^ ]+)/
        if (fmt_matcher) {
            out_fmt = fmt_matcher[0][1]
        }
    }
    if (out_fmt == "cram") {
        """
        pbmm2 \\
            align \\
            ${args} \\
            -j ${task.cpus} \\
            ${reference} \\
            ${reads} \\
            | samtools view \\
                -C \\
                -T ${fasta} \\
                -@ ${task.cpus} \\
                -o ${prefix}.cram -

        samtools index ${prefix}.cram
        """
    }
    else {
        """
        pbmm2 \\
            align \\
            ${args} \\
            -j ${task.cpus} \\
            ${reference} \\
            ${reads} \\
            ${prefix}.bam
        """
    }

    stub:
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_fmt = "bam"
    if (args2.contains('--output-fmt')) {
        def fmt_matcher = args2 =~ /--output-fmt\s+([^ ]+)/
        if (fmt_matcher) {
            out_fmt = fmt_matcher[0][1]
        }
    }
    if (out_fmt == "cram") {
        """
        touch ${prefix}.cram
        touch ${prefix}.cram.crai
        """
    }
    else {
        """
        touch ${prefix}.bam
        touch ${prefix}.bam.bai
        """
    }

    output:
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.cram"), optional: true, emit: cram
    tuple val(meta), path("*.bam.bai"), optional: true, emit: index
    tuple val(meta), path("*.cram.crai"), optional: true, emit: index_cram
    tuple val("${task.process}"), val("pbmm2"), eval("pbmm2 --version"), topic: versions, emit: versions
    tuple val("${task.process}"), val("samtools"), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
}
