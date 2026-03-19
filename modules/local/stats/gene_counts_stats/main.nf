process GENE_COUNTS_STATS {

    container 'docker.io/biocorecrg/multiqc:1.28'
    label 'process_low'

    input:
    path desc_file
    // desc.txt file 
    path annotated_file
    // norm_counts.genes file output by multiq_pca as emit: norm_counts
    val experiment
    // "rnaseq" or "smallrnaseq"
    val rna_type

    script:
    def rna_type_flag = rna_type != "" ? "-r ${rna_type}" : ""

    """
    RNA_summary.py --desc ${desc_file} --annotation ${annotated_file} --experiment ${experiment} ${rna_type_flag}
    """

    output:
    path ("*.csv"), emit: csv
}
