


process Fusioncatcher {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/fusion", mode: "${params.publishDirMode}"
 
    input:
    tuple val(meta),path(trim),path(db)
        
    output:
    tuple val(meta),path("${meta.lib}.fusion-catcher.txt"),  emit : fc_output
    tuple val(meta),path("${meta.lib}.summary_candidate_fusions.txt"), emit : fc_summary
    
    stub:
    """
    touch "${meta.lib}.fusion-catcher.txt"
    touch "${meta.lib}.summary_candidate_fusions.txt"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    # if running on biowulf SLURM
    if [ -d "/lscratch/${SLURM_JOB_ID}" ];then
        TMPDIR="/lscratch/${SLURM_JOB_ID}/${prefix}_STAR"
        if [ -d \${TMPDIR} ];then rm -rf \${TMPDIR};fi

        fusioncatcher.py \
            -p ${task.cpus} \
            -d ${db} \
            -i ${trim[0]},${trim[1]} \
            -o \${TMPDIR}
        
        cp \${TMPDIR}/final-list_candidate-fusion-genes.hg19.txt ${prefix}.fusion-catcher.txt
        cp \${TMPDIR}/summary_candidate_fusions.txt ${prefix}.summary_candidate_fusions.txt
    else
        mkdir fusioncatcher/
        fusioncatcher.py \
            -p ${task.cpus} \
            -d ${db} \
            -i ${trim[0]},${trim[1]} \
            -o fusioncatcher/

        cp fusioncatcher/final-list_candidate-fusion-genes.hg19.txt ${prefix}.fusion-catcher.txt
        cp fusioncatcher/summary_candidate_fusions.txt ${prefix}.summary_candidate_fusions.txt       
    fi

    """

}
