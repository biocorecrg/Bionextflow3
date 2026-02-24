process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/66dc96eff11ab80dfd5c044e9b3425f52d818847b9c074794cf0c02bfa781661/data' :
        'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(fasta)
    val bam_format
    val bam_index_extension
    val cigar_paf_format
    val cigar_bam

    output:
    tuple val(meta), path("*.paf")                       , optional: true, emit: paf
    tuple val(meta), path("*.bam")                       , optional: true, emit: bam
    tuple val(meta), path("*.cram")                       , optional: true, emit: cram
    tuple val(meta), path("*.bam.${bam_index_extension}"), optional: true, emit: index
    tuple val("${task.process}"), val("minimap2"), eval("minimap2 --version"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_fmt = "bam"
    if (args2.contains('--output-fmt')) {
        def fmt_matcher = args2 =~ /--output-fmt\s+([^ ]+)/
        if (fmt_matcher) out_fmt = fmt_matcher[0][1]
    }
    def samtools_ref = (out_fmt == "cram" && fasta) ? "--reference ${fasta}" : ""
    def bam_index = bam_index_extension ? "${prefix}.bam##idx##${prefix}.bam.${bam_index_extension} --write-index" : "${prefix}.${out_fmt}"
    def bam_output = bam_format ? "-a | samtools sort -@ ${task.cpus-1} -o ${bam_index} ${samtools_ref} ${args2}" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def samtools_reset_fastq = bam_input ? "samtools reset --threads ${task.cpus-1} $args3 $reads | samtools fastq --threads ${task.cpus-1} $args4 |" : ''
    def query = bam_input ? "-" : reads
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)

    """
    $samtools_reset_fastq \\
    minimap2 \\
        $args \\
        -t $task.cpus \\
        $target \\
        $query \\
        $cigar_paf \\
        $set_cigar_bam \\
        $bam_output


    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args2 = task.ext.args2 ?: ''
    def out_fmt = "bam"
    if (args2.contains('--output-fmt')) {
        def fmt_matcher = args2 =~ /--output-fmt\s+([^ ]+)/
        if (fmt_matcher) out_fmt = fmt_matcher[0][1]
    }
    def output_file = bam_format ? "${prefix}.${out_fmt}" : "${prefix}.paf"
    def bam_index = bam_index_extension ? "touch ${prefix}.bam.${bam_index_extension}" : ""
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)

    """
    touch $output_file
    ${bam_index}

    """
}
