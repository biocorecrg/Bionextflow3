process NANOCONSENSUS {
    tag "${meta.id}---${chrName}"
    label 'process_short_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker.io/biocorecrg/mop_consensus:0.1'
        : 'docker.io/biocorecrg/mop_consensus:0.1'}"

    input:

    tuple val(meta), path(Epi_Sample), path(Epi_IVT), path(baseQ), path(nanorms_si), path(nanorms_dt), path(nanorms_sd), val(chrName), val(chrStart), val(chrEnd)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(mod_annotation)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
	NanoConsensus.R -Epi_Sample ${Epi_Sample} \
	-Epi_IVT ${Epi_IVT} \
	-BaseQ ${baseQ} \
	-nanoRMS_SI ${nanorms_si} -nanoRMS_DT ${nanorms_dt} -nanoRMS_SD ${nanorms_sd}\
	-chr ${chrName} \
	-ini_pos ${chrStart} -fin_pos ${chrEnd} \
	-output ${prefix} \
	-fasta ${reference} --bed ${mod_annotation} ${args}
    """

    output:
    tuple val(meta), path("*-NanoConsensus_Scores.pdf"), emit: plots
    tuple val(meta), path("*_Supported_kmers.bedrmod"), emit: results
    tuple val(meta), path("BedRmod_tracks/*"), emit: bedrmods_tracks
    tuple val("${task.process}"), val("nanoconsensus"), eval("echo 2.0"), topic: versions
}
