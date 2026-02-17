process NANOSTAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanostat:1.6.0--pyhdfd78af_0' :
        'biocontainers/nanostat:1.6.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_file), path(index)

    output:
    tuple val(meta), path("*.txt")                 , emit: txt
    tuple val("${task.process}"), val("nanostat"), eval("NanoStat -v | cut -d ' ' -f 2"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
   
     script:
    def format_flag = input_file.name.endsWith('.cram') ? '--cram' : '--bam'

    """
    NanoStat \\
        $args \\
        -t $task.cpus \\
        ${format_flag} ${input_file} \\
        --tsv -o ./ \\
        -n ${meta.id}_nanostat.txt
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



    """
}
