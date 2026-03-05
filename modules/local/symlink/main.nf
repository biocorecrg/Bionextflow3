process SYMLINK {
    tag "$meta.id"
    label 'process_low'
    container 'docker.io/biocorecrg/almalinux-perlbrew-pyenv3'

    input:
    tuple val(meta), path(file, stageAs: 'input/*')
    
    output:
    tuple val(meta), path("${meta.id}*"), emit: linkedfile

    script:
    def extension = file.getExtension()
    """
    ln -s ${file} ${meta.id}.${extension}
    """

}
