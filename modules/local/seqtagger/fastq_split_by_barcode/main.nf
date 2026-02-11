process FASTQ_SPLIT_BY_BARCODE {
    tag "$meta.id"
    label 'gpu'

    container "docker://lpryszcz/seqtagger:1.0d"

    input:
	tuple val(meta), path(demux), path(fastq) 	
	
	output:
	tuple val(meta), path ("*.fq.gz"), emit: demux_fastq


	script:

    args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

	"""
		fastq_split_by_barcode.py -b 50 -i ${demux} -f ${fastq} -o ${prefix}
    """
}
