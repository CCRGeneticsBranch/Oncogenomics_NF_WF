

workflow HLA_calls {
    take: trimmed_fq
    main:
        HLAminer(trimmed_fq)
        Seq2HLA(trimmed_fq)
        MergeHLA(Seq2HLA.out.combine(HLAminer.out, by:0))
    emit:
        MergeHLA.out
}


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

process Seq2HLA {
    tag { dataset_id }

    publishDir "$params.resultsdir/$dataset_id/HLA", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2)
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}-ClassI.HLAgenotype4digits")

    stub:
    """
    touch "${dataset_id}-ClassI.HLAgenotype4digits"
    """

    shell:
    '''
    set -exo pipefail

    seq2HLA \
        references/ \
        -1 !{r1} \
        -2 !{r2} \
        -p !{task.cpus} \
        -r "!{dataset_id}"
    '''
}


process MergeHLA {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/HLA", mode: "${params.publishDirMode}"
    input:
    tuple val(dataset_id),
        path(seq2hla),
        path(hlaminer)

    output:
	path "${dataset_id}.Calls.txt"

    stub:
    """
    touch "${dataset_id}.Calls.txt"
    """


    shell:
    '''
    set -exo pipefail

    consensusHLA.pl !{hlaminer} !{seq2hla} |sort > !{dataset_id}.Calls.txt

    '''
}


