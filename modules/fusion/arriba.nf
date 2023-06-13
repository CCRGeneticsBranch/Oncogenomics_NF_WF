process Arriba {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/fusion", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(trim),path(reffa),path(star_genomeIndex),path(gtf)
 
    output:
    tuple val(meta),path("${meta.lib}.arriba-fusion.txt")
 
    stub:
    """
      touch "${meta.lib}.arriba-fusion.txt"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    set -exo pipefail
    
    # if running on biowulf SLURM
    if [ -d "/lscratch/${SLURM_JOB_ID}" ];then
        TMPDIR="/lscratch/${SLURM_JOB_ID}/${prefix}_STAR"
        if [ -d \${TMPDIR} ];then rm -rf \${TMPDIR};fi
            
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
            --outTmpDir \${TMPDIR} \
            --outFileNamePrefix ${prefix}.arriba.
        
        samtools sort -@ ${task.cpus} -T \$TMPDIR -o ${prefix}.arriba.Aligned.sortedByCoords.out.bam -O BAM ${prefix}.arriba.Aligned.out.bam
        
    else
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
        
        samtools sort -@ ${task.cpus} -o ${prefix}.arriba.Aligned.sortedByCoords.out.bam -O BAM ${prefix}.arriba.Aligned.out.bam
    fi
    
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

    mv ${prefix}.fusions.tsv ${prefix}.arriba-fusion.txt
    """

}
