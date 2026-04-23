process METHODS_SECTION {

    label 'process_low'
    container 'ghcr.io/biocorecrg/bionextflow3/biocorecrg_methods:edge'

    input:
    path pipeline_paths, stageAs: "?/*"
    path(params_file)
    path(template)
    path(citations_file)
    path(nf_core)

    output:
  //  path("*_output.yml"), emit: methods_section
    path("methods_description_mqc.yml"), emit: methods_mqc

    script:
    // Check if optional file exists and build the flag conditionally
    def citations_arg = citations_file ? "-c ${citations_file}" : ''

    """
    get_doi_from_meta.py -p ./
    addCitationFromYaml.js --input dois.yml --output methods_data.yml 
    pipeline_methods.py -d methods_data.yml -m ${params_file} -t template.yml ${citations_arg}
    """
}
  