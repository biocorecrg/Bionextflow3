process EPINANO_PLOTS {
    tag "$meta.id"
    label 'process_middle'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/biocorecrg/epinanoplots:0.1' :
        'docker.io/biocorecrg/epinanoplots:0.1' }"

    input:
    tuple val(meta), path(alnfile), path(alnindex)
    val(mode)

   // input:
   // tuple val(sampleIDA), val(sampleIDB), path(per_site_varA), path(per_site_varB)


    output:
    tuple val(meta), path("*.pdf"), emit: pdf
  
    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    script:
	"""
	#epinano_scatterplot.R ${per_site_varA} ${sampleIDA} ${per_site_varB} ${sampleIDB} ${mode}
    """


}
