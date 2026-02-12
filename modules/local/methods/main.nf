process METHODS_SECTION {

    label 'process_low'
    container 'biocorecrg/methods:0.1'

    input:
    path pipeline_paths, stageAs: "?/*"
    path(params_file)
    path(template)
    path(nf_core)

    output:
  //  path("*_output.yml"), emit: methods_section
    path("methods_description_mqc.yml"), emit: methods_mqc

    script:
    """
    get_doi_from_meta.py -p ./
    addCitationFromYaml.js --input dois.yml --output methods_data.yml 
    pipeline_methods.py -d methods_data.yml -m methods.yaml -t template.yml
 """
}
