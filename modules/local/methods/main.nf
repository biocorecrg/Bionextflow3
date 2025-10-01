process METHODS_SECTION {

    label 'process_low'

    input:
    path(pipeline_path)
    path(params_file)

    output:
    path("*_output.yml"), emit: methods_section
    path("methods_description_mqc.yml"), emit: methods_mqc

    script:
    """
    pipeline_methods.py -p ${pipeline_path}
    addCitationFromYaml.js --input methods_output_pre.yml --output methods_output.yml
    smallrnaseq_methods.py --methods_output methods_output.yml --params ${params_file}
    """
}
