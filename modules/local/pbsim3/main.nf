process PBSIM3 {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/biocontainers/pbsim3:3.0.5--h9948957_2'
        : 'biocontainers/pbsim3:3.0.5--h9948957_2'}"

    input:
    tuple val(meta), path(genome)
    tuple val(meta2), path(method_file)
    val sim_method
    val depth

    output:
    tuple val(meta), path("*.fastq"), emit: fastq
    tuple val("${task.process}"), val("pbsim3"), eval("pbsim --version | head -n 1"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // sim_method is the pbsim3 method name directly: qshmm, errhmm, or sample
    def valid_methods = ['qshmm', 'errhmm', 'sample']
    if (!valid_methods.contains(sim_method)) {
        error "Unknown sim_method '${sim_method}'. Must be one of: ${valid_methods.join(', ')}"
    }
    def method_arg = "--method ${sim_method} --${sim_method} ${method_file}"

    def fasta_name = genome.name.endsWith('.gz') ? genome.baseName : genome.name
    def decompress_cmd = genome.name.endsWith('.gz') ? "gunzip -k ${genome}" : ''

    """
    ${decompress_cmd}
    pbsim \\
        --strategy wgs \\
        ${method_arg} \\
        --depth ${depth} \\
        --genome ${fasta_name} \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_0001.fastq
    """
}
