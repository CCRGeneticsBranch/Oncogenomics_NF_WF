process Strandedness {
     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

     input:
        tuple val(meta),path(G_bam),path(G_bai),path(gtf_sorted),path(gtf_index)
 
     output:
     tuple val(meta),path("${meta.lib}_strandedness.txt")
 
     script:
      def args = task.ext.args   ?: ''
      def prefix   = task.ext.prefix ?: "${meta.lib}"

     """
     ngsderive strandedness -g $gtf_sorted $G_bam -n 10000 > ${prefix}_strandedness.txt

     """
}

process Rsem {
        tag "$meta.lib"

        publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/RSEM_ENS", mode: 'copy',pattern: "${meta.lib}*"

        input:
        tuple val(meta), path(T_bam), path(strandedness), path(genomeIndex)
 
        output:
        tuple val(meta), path("${meta.lib}.genes.results"), emit: rsem_genes
        tuple val(meta), path("${meta.lib}.isoforms.results"), emit: rsem_isoforms
        path "versions.yml"             , emit: versions

        script:
        def args = task.ext.args   ?: ''
        def prefix   = task.ext.prefix ?: "${meta.lib}"

        """
        STRAND=`strandedness.py ${prefix}_strandedness.txt rsem`
        echo "strandedness is  \$STRAND"
        rsem-calculate-expression --no-bam-output --paired-end --strandedness \$STRAND -p ${task.cpus} --estimate-rspd --bam $T_bam ${genomeIndex}/rsem_1.3.2 ${prefix}
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
                RSEM: \$(rsem-calculate-expression --version |cut -f4 -d " ")
        END_VERSIONS
        """

}

