process POD5_FILTER {
    tag "$meta.id"
    label 'process_medium'

    container {
        def arch = System.getProperty("os.arch")
        if (arch == "aarch64") {
            "community.wave.seqera.io/library/pod5:0.3.33--05f974e602a2b794"  // Apple Silicon
        } else {
            "community.wave.seqera.io/library/pod5:0.3.33--762a6be5bc773897"  // Intel x86-64
        }
    }
    input:
    tuple val(meta), path(list), path(pod5)
    
    output:
    tuple val(meta), path("*.pod5")          , emit: pod5
    tuple val("${task.process}"), val("pod5"), eval("pod5 --version 2>&1"), topic: versions, emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    pod5 filter -t ${task.cpus}  -M -D -i ${list} -o ${meta.id}.pod5 ${pod5}
    """
}
