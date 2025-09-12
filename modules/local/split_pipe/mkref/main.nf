process SPLITPIPE_MKREF {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'biocorecrg/spipe:1.6.1'

    input:
    tuple val(meta),   path(genome)
    path(annotation)
    
    output:
    tuple val(meta), path("${meta.id}")                       , emit: index
    path  "versions.yml"                                      , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split-pipe: \$(split-pipe --version | cut -d ' ' -f 2 | sed 's/v//g')
END_VERSIONS
    """

}
