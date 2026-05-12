process EVENTUAL_GZIP_FASTQ {
    label "process_medium"
    tag "${meta.id}"

    container 'quay.io/biocontainers/gzip:1.11'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(outputs), emit: zipped
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | head -n 1 | sed "s/^.*gzip //; s/ .*\$//"'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (meta.single_end) {
        outputs = reads.extension == 'gz' ? reads.name : "${reads.name}.gz"
        if (reads.extension == 'gz') {
            """
            :
            """
        } else {
            """
            gzip -c ${reads} > ${outputs}
            """
        }
    } else {
        def r1 = reads[0]
        def r2 = reads[1]
        def out1 = r1.extension == 'gz' ? r1.name : "${r1.name}.gz"
        def out2 = r2.extension == 'gz' ? r2.name : "${r2.name}.gz"
        outputs = [out1, out2]
        def cmd1 = r1.extension == 'gz' ? ":" : "gzip -c ${r1} > ${out1}"
        def cmd2 = r2.extension == 'gz' ? ":" : "gzip -c ${r2} > ${out2}"
        """
        ${cmd1}
        ${cmd2}
        """
    }

    stub:
    if (meta.single_end) {
        outputs = reads.extension == 'gz' ? reads.name : "${reads.name}.gz"
    } else {
        def out1 = reads[0].extension == 'gz' ? reads[0].name : "${reads[0].name}.gz"
        def out2 = reads[1].extension == 'gz' ? reads[1].name : "${reads[1].name}.gz"
        outputs = [out1, out2]
    }
    """
    touch ${outputs instanceof List ? outputs.join(' ') : outputs}
    """
}
