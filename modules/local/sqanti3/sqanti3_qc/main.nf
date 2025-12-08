process SQANTI3_QC {
    tag "${meta.id}"
    label 'process_medium'

    container 'anaconesalab/sqanti3:v5.5.4'

    input:
    tuple val(meta), path(gtf)               // GTF file to evaluate
    tuple val(meta2), path(fasta), path(fai) // Reference genome FASTA and .fai
    tuple val(meta3), path(reference_gtf)    // Reference annotation GTF

    output:
    tuple val(meta), path("sqanti3_results/*_classification.txt")          , emit: classification, optional: true
    tuple val(meta), path("sqanti3_results/*_junctions.txt")               , emit: junctions, optional: true
    tuple val(meta), path("sqanti3_results/*_corrected.gtf")               , emit: corrected_gtf, optional: true
    tuple val(meta), path("sqanti3_results/*_corrected.fasta")             , emit: corrected_fasta, optional: true
    tuple val(meta), path("sqanti3_results/*.pdf")                         , emit: plots, optional: true
    tuple val(meta), path("sqanti3_results/*.html")                        , emit: report, optional: true
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # Run SQANTI3 qc:
    /conda/miniconda3/envs/sqanti3/bin/python /opt2/sqanti3/5.5.4/SQANTI3-5.5.4/sqanti3_qc.py \\
    --isoforms ${gtf} \\
    --refFasta ${fasta} \\
    --refGTF ${reference_gtf} \\
    --cpus ${task.cpus} \\
    ${args}
    
    # Generate versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_qc.py --version 2>&1 | grep -oP 'SQANTI3 v\\K[0-9.]+' || echo "5.2.2")
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """


}
