process Fusioncatcher{
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/fusioncatcher", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2),
        path(db)
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}.final-list_candidate-fusion-genes.hg19.txt"),
        path("${dataset_id}.summary_candidate_fusions.txt")

    stub:
    """
    touch "${dataset_id}.final-list_candidate-fusion-genes.hg19.txt"
    touch "${dataset_id}.summary_candidate_fusions.txt"
    """

    shell:
    '''
    # if running on biowulf SLURM
    if [ -d "/lscratch/${SLURM_JOB_ID}" ];then
        TMPDIR="/lscratch/${SLURM_JOB_ID}/!{dataset_id}_STAR"
        if [ -d ${TMPDIR} ];then rm -rf ${TMPDIR};fi

        fusioncatcher.py \
            -p !{task.cpus} \
            -d !{db} \
            -i !{r1},!{r2} \
            -o ${TMPDIR}
        
        cp ${TMPDIR}/final-list_candidate-fusion-genes.hg19.txt !{dataset_id}.final-list_candidate-fusion-genes.hg19.txt
        cp ${TMPDIR}/summary_candidate_fusions.txt !{dataset_id}.summary_candidate_fusions.txt
    else
        mkdir fusioncatcher/
        fusioncatcher.py \
            -p !{task.cpus} \
            -d !{db} \
            -i !{r1},!{r2}
            -o fusioncatcher/
        
        cp fusioncatcher/final-list_candidate-fusion-genes.hg19.txt !{dataset_id}.final-list_candidate-fusion-genes.hg19.txt
        cp fusioncatcher/summary_candidate_fusions.txt !{dataset_id}.summary_candidate_fusions.txt
    fi

    '''

}
