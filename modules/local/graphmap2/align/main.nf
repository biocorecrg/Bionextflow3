process GRAPHMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0' :
        'biocontainers/graphmap:0.6.3--he513fc3_0' }"

    input:
    tuple val(meta), path(reads)
    path  fasta
    path  index

    output:
    tuple val(meta), path("*.cram"), emit: cram
    tuple val("${task.process}"), val('graphmap2'), eval("graphmap2 align 2>&1 | head -1 | sed 's/.*Version: v//; s/ .*//' || true"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    graphmap2 \\
        align \\
        -t $task.cpus \\
        -r $fasta \\
        -i $index \\
        -d $reads \\
        -o ${prefix}.sam \\
        $args

    samtools view -@ ${task.cpus} -F4 -T $fasta -hSC ${prefix}.sam > ${prefix}.cram

    rm ${prefix}.sam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.cram
    """

}
