process HLAminer {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/HLAminer", mode: "${params.publishDirMode}"
    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2)
    
    output:
        path "HLAminer_HPTASR.csv"

    stub:
    """
    touch "HLAminer_HPTASR.csv"
    """


    shell:
    '''
    set -exo pipefail

    
    echo !{r1} >patient.fof
    echo !{r2} >>patient.fof


    HPTASRrnaseq_classI.sh .

    '''
}
