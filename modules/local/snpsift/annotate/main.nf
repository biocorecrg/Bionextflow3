process SNPSIFT_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3d8ec79a01bcc86a5ce258c66fc18e48c1826aebc7e7114454757919162ff9e6/data'
        : 'community.wave.seqera.io/library/snpsift:5.4.0c--6546f37f72acfb46'}"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    path databases
    path databases_idx

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_list = databases instanceof List ? databases : [databases]
    def commands = []

    db_list.eachWithIndex { db, idx ->
        // Extract database base name and clean it up (e.g. clinvar.vcf.gz -> CLINVAR)
        def db_name = db.name.replaceAll(/(\.vcf)?(\.gz)?(\.tbi)?$/, '')
        def db_prefix = "${db_name.toUpperCase().replaceAll(/[^A-Z0-9_]/, '_')}_"
        def current_args = "${args} -name ${db_prefix}"

        if (idx == 0) {
            commands << "SnpSift annotate ${current_args} ${db} ${vcf}"
        }
        else {
            commands << "SnpSift annotate ${current_args} ${db}"
        }
    }
    def final_command = commands.join(" \\\n        | ") + " > ${prefix}.vcf"
    """
    mkdir -p \$PWD/tmp
    export JAVA_OPTS="-Djava.io.tmpdir=\$PWD/tmp"
    ${final_command}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('snpsift'), eval("SnpSift split -h 2>&1 | sed 's/^.*version //' | sed 's/(.*//' | sed 's/ .*//'"), topic: versions, emit: versions_snpsift
}
