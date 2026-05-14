process MODELS_TO_MULTIQC {
    tag "${meta.id}"
    label 'small'

    container "docker.io/nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"

    input:
    tuple val(meta), path(model_folder)

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''


    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "# id: 'dorado_models_info'" > dorado_info_mqc.tsv
    echo "# section_name: 'Dorado Models Information'" >> dorado_info_mqc.tsv
    echo "# plot_type: 'table'" >> dorado_info_mqc.tsv
    echo -e "Number\tModel used" >> dorado_info_mqc.tsv
    ls ${model_folder} | awk '{num++; print num"\t"\$0}' >> dorado_info_mqc.tsv
    """

    output:
    tuple val(meta), path("dorado_info_mqc.tsv"), emit: mqc
}
