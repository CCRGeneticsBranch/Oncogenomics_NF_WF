


process Fusioncatcher{
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/fusion", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        val(library),
        path(r1), 
        path(r2),
        path(db)
    
    output:
    tuple val("${dataset_id}"),
        val("${library}"),
        path("${library}.fusion-catcher.txt"),
        path("${library}.summary_candidate_fusions.txt")

    stub:
    """
    touch "${library}.fusion-catcher.txt"
    touch "${library}.summary_candidate_fusions.txt"
    """

    shell:
    '''
    # if running on biowulf SLURM
    if [ -d "/lscratch/${SLURM_JOB_ID}" ];then
        TMPDIR="/lscratch/${SLURM_JOB_ID}/!{library}_STAR"
        if [ -d ${TMPDIR} ];then rm -rf ${TMPDIR};fi

        fusioncatcher.py \
            -p !{task.cpus} \
            -d !{db} \
            -i !{r1},!{r2} \
            -o ${TMPDIR}
        
        cp ${TMPDIR}/final-list_candidate-fusion-genes.hg19.txt !{library}.fusion-catcher.txt
        cp ${TMPDIR}/summary_candidate_fusions.txt !{library}.summary_candidate_fusions.txt
    else
        mkdir fusioncatcher/
        fusioncatcher.py \
            -p !{task.cpus} \
            -d !{db} \
            -i !{r1},!{r2} \
            -o fusioncatcher/

        cp fusioncatcher/final-list_candidate-fusion-genes.hg19.txt !{library}.fusion-catcher.txt
        cp fusioncatcher/summary_candidate_fusions.txt !{library}.summary_candidate_fusions.txt       
    fi

    '''

}
