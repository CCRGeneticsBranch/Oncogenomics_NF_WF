process Strandedness {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/qc", mode: "${params.publishDirMode}"

     input:
        tuple val(dataset_id),
        val(library),
        path(T_bam),
        path(G_bam),
        path(G_bai),
        path(chimeric_junction),
        path(gtf_sorted),
        path(gtf_index)

     output:
     tuple val("${dataset_id}"),
       val("${library}"),
       path("${library}_strandedness.txt")

     script:
     """
     ngsderive strandedness -g $gtf_sorted $G_bam -n 10000 > ${library}_strandedness.txt

     """
}

process Rsem {
        tag { dataset_id }

        publishDir "$params.resultsdir/$dataset_id/${params.casename}/${library}/RSEM_ENS", mode: 'copy'

        input:

        tuple val(dataset_id),
        val(library),
        path(T_bam),
        path(G_bam),
        path(G_bai),
        path(chimeric_junction),
        path(strandedness),
        path(genomeIndex)

        output:
        tuple val("${dataset_id}"),
                val(library),
                path("${library}.genes.results")


        script:
        """
        STRAND=`strandedness.py ${library}_strandedness.txt rsem`
        echo "strandedness is  \$STRAND"
        rsem-calculate-expression --no-bam-output --paired-end --strandedness \$STRAND -p ${task.cpus} --estimate-rspd --bam $T_bam ${genomeIndex}/rsem_1.3.2 $library
        """

}

