process DORADO_DEMUX {
    tag "$meta.id"
    label 'process_high'

    container "ontresearch/dorado:sha00aa724a69ddc5f47d82bd413039f912fdaf4e77"

    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.bam")          , emit: basecalled_bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    def args = task.ext.args ?: ''
 
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dorado demux --emit-summary $args --threads ${task.cpus} --output-dir ./ ${bam}
    gzip *_summary.txt
    # Avoid the folder structure from dorado
    mv */*/*/bam_pass/*/*.bam .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -n1 | cut -d "-" -f 2 | cut -d "+" -f 1) 
END_VERSIONS
    """
}
