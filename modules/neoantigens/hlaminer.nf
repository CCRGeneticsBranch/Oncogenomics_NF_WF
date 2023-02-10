process HLAminer {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/HLA", mode: "${params.publishDirMode}"
    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2)
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}_HLAminer_HPTASR.csv")


    stub:
    """
    touch "${dataset_id}_HLAminer_HPTASR.csv"
    """


    shell:
    '''
    set -exo pipefail

    
    echo !{r1} >patient.fof
    echo !{r2} >>patient.fof


    HPTASRrnaseq_classI.sh .
    mv HLAminer_HPTASR.csv !{dataset_id}_HLAminer_HPTASR.csv
    '''
}
