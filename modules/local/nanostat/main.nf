process NANOSTAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanostat:1.6.0--pyhdfd78af_0' :
        'biocontainers/nanostat:1.6.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam_file), path(bai_file)

    output:
    tuple val(meta), path("*.txt")                 , emit: txt
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
   
    """
    NanoStat \\
        $args \\
        -t $task.cpus \\
        --bam ${bam_file} \\
        --tsv -o ./ \\
        -n ${meta.id}_nanostat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanostat: \$(echo \$(NanoStat -v 2>&1) | sed 's/^.*NanoStat //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch LengthvsQualityScatterPlot_dot.html
    touch LengthvsQualityScatterPlot_kde.html
    touch NanoPlot-report.html
    touch NanoStats.txt
    touch Non_weightedHistogramReadlength.html
    touch Non_weightedLogTransformed_HistogramReadlength.html
    touch WeightedHistogramReadlength.html
    touch WeightedLogTransformed_HistogramReadlength.html
    touch Yield_By_Length.html


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
