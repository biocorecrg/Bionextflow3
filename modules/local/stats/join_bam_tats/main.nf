process JOIN_BAM_STATS {

    tag "joining aln stats"
    label     'low'
    shell '/bin/bash'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1a35167f7a491c7086c13835aaa74b39f1f43979:9254eac8981f615fb6c417fa44e77c3b44bc3abd-0' :
        'quay.io/biocontainers/mulled-v2-1a35167f7a491c7086c13835aaa74b39f1f43979:a7b00ff483a30f0a985d9e0d4da1f5762af68cd6-0' }"

    input:
    file "alnqc_*"

    output:
    path("alnQC_mqc.txt"), emit: join_stats

    
    """
    echo '# id: alnQC
# plot_type: \'table\'
# section_name: \'Alignment QC\' ' > alnQC_mqc.txt
    cat alnqc_* | head -n 1| sed s@#@@g >> alnQC_mqc.txt
    cat alnqc_* | grep -v "#" >> alnQC_mqc.txt
     """
}
