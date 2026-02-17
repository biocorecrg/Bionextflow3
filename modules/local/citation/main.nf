process CITATION_JS {
    tag "${meta.id}"
    label 'process_low'

    // Use container or otherwise use wave using Dockerfile available
    container "docker.io/biocorecrg/citation-js:0.7.18"

    input:
    tuple val(meta), val(doi), val(format), val(citation_template)

    output:
    tuple val(meta), val(output)                    , emit: output
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
    retrieveFromDOI.js --format ${format} --template ${citation_template} ${doi}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        citation-js: 0.7.18
    END_VERSIONS
    """

}
