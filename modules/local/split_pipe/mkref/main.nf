process SPLITPIPE_MKREF {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/biocorecrg/spipe:1.6.1'

    input:
    tuple val(meta),   path(genome)
    path(annotation)
    
    output:
    tuple val(meta), path("${meta.id}")                       , emit: index
    tuple val("${task.process}"), val('split-pipe'), eval('split-pipe --version | cut -d " " -f 2 | sed "s/v//g"'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
    split-pipe \
	   --mode mkref ${args} \
	   --genome_name ${meta.id}  \
	   --fasta ${genome} \
       --genes ${annotation} \
       --nthreads ${task.cpus} \
	   --output_dir ./${meta.id}

	# Fix for missing version
	sed -i s/'"."'/'"1.6.1"'/ ${meta.id}/process/mkref_def.json

    """

}
