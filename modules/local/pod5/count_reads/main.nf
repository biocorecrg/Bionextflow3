process POD5_COUNT_READS {
    tag "$meta.id"
    label 'process_low'

    container {
        def arch = System.getProperty("os.arch")
        if (arch == "aarch64") {
            "community.wave.seqera.io/library/pod5:0.3.33--05f974e602a2b794"
        } else {
            "community.wave.seqera.io/library/pod5:0.3.33--762a6be5bc773897"
        }
    }

    input:
    tuple val(meta), path(pod5)

    output:
    tuple val(meta), path(pod5), path("reads.txt"), emit: count
    tuple val("${task.process}"), val("pod5"), eval("pod5 --version 2>&1"), topic: versions, emit: versions

    script:
    """
    python3 -c "import pod5 as p5; print(p5.Reader('${pod5}').num_reads)" > reads.txt
    """
}
