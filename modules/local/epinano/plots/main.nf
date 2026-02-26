process EPINANO_PLOTS {
    tag "$meta.id vs $meta2.id"
    label 'process_middle'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/biocorecrg/epinanoplots:0.1' :
        'docker.io/biocorecrg/epinanoplots:0.1' }"

    input:
    tuple val(meta), path(per_site_varA)
    tuple val(meta2), path(per_site_varB)
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
	epinano_scatterplot.R ${per_site_varA} ${meta.id} ${per_site_varB} ${meta2.id} ${mode}
    """


}
