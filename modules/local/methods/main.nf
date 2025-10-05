process METHODS_SECTION {

    label 'process_low'
    container 'biocorecrg/methods:0.1'

    input:
    path(pipeline_paths)
    path(params_file)
    path(template)
    path(nf_core)

    output:
  //  path("*_output.yml"), emit: methods_section
    path("methods_description_mqc.yml"), emit: methods_mqc

    script:
    """
    pipeline_methods.py -p ${projectDir}
    addCitationFromYaml.js --input methods_output_pre.yml --output methods_output.yml
   smallrnaseq_methods.py -v methods_output.yml -m methods.yaml -t template.yml """
}
