process Arriba{
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/fusion", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        val(library),
        path(r1), 
        path(r2),
        path(reffa),
        path(star_genomeIndex),
        path(gtf)

    output:
    tuple val("${dataset_id}"),
        val("${library}"),
        path("${library}.arriba-fusion.txt")

    stub:
    """
      touch "${library}.arriba-fusion.txt"
    """

    shell:
    '''
    set -exo pipefail
    
    # if running on biowulf SLURM
    if [ -d "/lscratch/${SLURM_JOB_ID}" ];then
        TMPDIR="/lscratch/${SLURM_JOB_ID}/!{library}_STAR"
        if [ -d ${TMPDIR} ];then rm -rf ${TMPDIR};fi
            
        STAR --genomeDir !{star_genomeIndex} \
            --readFilesIn !{r1} !{r2} \
            --readFilesCommand zcat \
            --runThreadN !{task.cpus} \
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
            --outTmpDir ${TMPDIR} \
            --outFileNamePrefix "!{library}.arriba."
        
        samtools sort -@ !{task.cpus} -T $TMPDIR -o !{library}.arriba.Aligned.sortedByCoords.out.bam -O BAM !{library}.arriba.Aligned.out.bam
        
    else
        STAR --genomeDir !{star_genomeIndex} \
            --readFilesIn !{r1} !{r2} \
            --readFilesCommand zcat \
            --runThreadN !{task.cpus} \
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
            --outFileNamePrefix "!{library}.arriba."
        
        samtools sort -@ !{task.cpus} -o !{library}.arriba.Aligned.sortedByCoords.out.bam -O BAM !{library}.arriba.Aligned.out.bam
    fi
    
    # Run arriba
    arriba \
        -x !{library}.arriba.Aligned.out.bam \
        -o !{library}.fusions.tsv \
        -O !{library}.fusions.discarded.tsv \
        -a !{reffa} \
        -g !{gtf} \
        -b /opt2/arriba_v2.3.0/database/blacklist_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz \
        -p /opt2/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3
        
    # index file
    samtools index -@ !{task.cpus} !{library}.arriba.Aligned.sortedByCoords.out.bam
    
    # process fusions    
    nfusions="$(wc -l !{library}.fusions.tsv | awk '{print \$1}')"
        
    if [ "$nfusions" -gt "1" ];then
        draw_fusions.R \
            --fusions=!{library}.fusions.tsv  \
            --alignments=!{library}.arriba.Aligned.sortedByCoords.out.bam \
            --output=!{library}.fusions.pdf \
            --annotation=!{gtf} \
            --cytobands=/opt2/arriba_v2.3.0/database/cytobands_hg19_hs37d5_GRCh37_v2.3.0.tsv \
            --proteinDomains=/opt2/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3
    else
        touch !{library}.fusions.pdf
    fi

    mv !{library}.fusions.tsv !{library}.arriba-fusion.txt
    '''

}
