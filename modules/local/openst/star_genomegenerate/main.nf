process STAR_GENOMEGENERATE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/biocorecrg/spacemake:0.9.1' :
        'docker.io/biocorecrg/spacemake:0.9.1' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("star_index")  , emit: species_folder
    tuple val("${task.process}"), val('STAR'), eval('/opt/conda/bin/STAR-avx2 --version'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"


    """
 	/opt/conda/bin/STAR-avx2 --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir star_index --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch species_data

    """
}
