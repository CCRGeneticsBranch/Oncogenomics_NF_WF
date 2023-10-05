process Arriba {
    tag "$meta.lib"
    scratch true
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/fusion", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(trim),path(reffa),path(star_genomeIndex),path(gtf)
 
    output:
    tuple val(meta),path("${meta.lib}.arriba-fusion.txt"), emit: arriba_fusion
    tuple val(meta),path("arriba_out/${meta.lib}.fusions.discarded.tsv"), emit: arriba_discarded
    tuple val(meta),path("arriba_out/${meta.lib}.fusions.pdf") , emit: arriba_pdf
 
    stub:
    """
      touch "${meta.lib}.arriba-fusion.txt"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    set -exo pipefail
    
    TMP=tmp/
    mkdir \$TMP
    trap 'rm -rf "\$TMP"' EXIT
            
        STAR --genomeDir ${star_genomeIndex} \
            --readFilesIn ${trim[0]} ${trim[1]} \
            --readFilesCommand zcat \
            --runThreadN ${task.cpus} \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outBAMcompression 0 \
            --outFilterMultimapNmax 50 \
            --peOverlapNbasesMin 10 \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimSegmentMin 10 \
            --chimOutType WithinBAM HardClip \
            --chimJunctionOverhangMin 10 \
            --chimScoreDropMax 30 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50 \
            --outFileNamePrefix ${prefix}.arriba.
        
        samtools sort -@ ${task.cpus} -T \$TMP -o ${prefix}.arriba.Aligned.sortedByCoords.out.bam -O BAM ${prefix}.arriba.Aligned.out.bam
        
    
    # Run arriba
    arriba \
        -x ${prefix}.arriba.Aligned.out.bam \
        -o ${prefix}.fusions.tsv \
        -O ${prefix}.fusions.discarded.tsv \
        -a ${reffa} \
        -g ${gtf} \
        -b /opt2/arriba_v2.3.0/database/blacklist_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz \
        -p /opt2/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3
        
    # index file
    samtools index -@ ${task.cpus} ${prefix}.arriba.Aligned.sortedByCoords.out.bam
    
    # process fusions    
    mkdir arriba_out
    NFUSIONS=`wc -l ${prefix}.fusions.tsv`
        
    if [ "\$NFUSIONS" -gt "1" ];then
        draw_fusions.R \
            --fusions=${prefix}.fusions.tsv  \
            --alignments=${prefix}.arriba.Aligned.sortedByCoords.out.bam \
            --output=${prefix}.fusions.pdf \
            --annotation=${gtf} \
            --cytobands=/opt2/arriba_v2.3.0/database/cytobands_hg19_hs37d5_GRCh37_v2.3.0.tsv \
            --proteinDomains=/opt2/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3
    else
        touch ${prefix}.fusions.pdf
    fi

    cp ${prefix}.fusions.tsv ${prefix}.arriba-fusion.txt
    
    mv ${prefix}.fusions.* ./arriba_out
    
    """

}
