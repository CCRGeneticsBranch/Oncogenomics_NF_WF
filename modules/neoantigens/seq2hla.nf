process Seq2HLA {
    tag { dataset_id }

    publishDir "$params.resultsdir/$dataset_id/Seq2HLA", mode: 'copy'

    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2),
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}.Aligned.toTranscriptome.out.bam"),
        path("${dataset_id}.Aligned.sortedByCoord.out.bam"),
        path("${dataset_id}.Aligned.sortedByCoord.out.bam.bai")

    script:
    """
    set -exo pipefail

    python seq2HLA.py \
        seq2HLA/ \
        -1 $r1 \
        -2 $r2 \
        -p ${task.cpus} \
        -r HLA/seq2HLA/{wildcards.sample}
    """
}