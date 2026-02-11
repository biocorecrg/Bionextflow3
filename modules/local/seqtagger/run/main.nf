process SEQTAGGER {
    tag "$meta.id"
    label 'gpu'

    container "docker://lpryszcz/seqtagger:1.0d"

    input:
    tuple val(meta), path(pod5)
    path  models
    
    output:
 	tuple val(meta), path("*_demux.tsv.gz"), emit: demux_files
	tuple val(meta), path("*.boxplot.pdf") , emit: demux_boxplot
    
    when:
    task.ext.when == null || task.ext.when

    script:
 
    args = task.ext.args ?: ''
 
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
        mkdir tmp
        export MPLCONFIGDIR=$PWD/tmp
    	run ${args} -r -i ./ -o temp_output -t ${task.cpus}
    	mv temp_output/..demux.tsv.gz ${prefix}_demux.tsv.gz
    	mv temp_output/..demux.tsv.gz.boxplot.pdf ${prefix}.boxplot.pdf
      
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtagger: \$(run --version 2>&1) 
END_VERSIONS
    """
}
