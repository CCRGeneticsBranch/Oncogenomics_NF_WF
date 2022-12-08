process Arriba{
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/Arriba", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2),
        path(reffa),
        path(genomeIndex),
        path(gtf)
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}.fusions.tsv"),
        path("${dataset_id}.fusions.discarded.tsv"),
        path("${dataset_id}.fusions.pdf")

    stub:
    """
    touch "${dataset_id}.fusions.tsv"
    touch "${dataset_id}.fusions.discarded.tsv"
    touch "${dataset_id}.fusions.pdf"
    """

    shell:
    '''
set -exo pipefail
if [ -d /lscratch/${SLURM_JOB_ID} ];then
    TMPDIR="/lscratch/${SLURM_JOB_ID}/!{dataset_id}_Arriba"
else
    TMPDIR="/dev/shm/!{dataset_id}_Arriba"
fi
if [ -d ${TMPDIR} ];then rm -rf ${TMPDIR};fi

STAR \
    --runThreadN !{task.cpus} \
    --genomeDir !{genomeIndex} \
    --readFilesIn !{r1} !{r2} \
    --readFilesCommand zcat \
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
    --outTmpDir $TMPDIR \
    --outFileNamePrefix "!{dataset_id}.arriba."

arriba \
    -x !{dataset_id}.arriba.Aligned.out.bam \
    -o !{dataset_id}.fusions.tsv \
    -O !{dataset_id}.fusions.discarded.tsv \
    -a !{reffa} \
    -g !{gtf} \
    -b /opt2/arriba_v2.3.0/database/blacklist_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz \
    -p /opt2/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3

    samtools sort -@ !{task.cpus} -T $TMPDIR -o !{dataset_id}.arriba.Aligned.sortedByCoords.out.bam -O BAM !{dataset_id}.arriba.Aligned.out.bam
    samtools index -@ !{task.cpus} !{dataset_id}.arriba.Aligned.sortedByCoords.out.bam

nfusions="$(wc -l !{dataset_id}.fusions.tsv | awk '{print \$1}')"

if [ "$nfusions" -gt "1" ];then
draw_fusions.R \
    --fusions=!{dataset_id}.fusions.tsv  \
    --alignments=!{dataset_id}.arriba.Aligned.sortedByCoords.out.bam \
    --output=!{dataset_id}.fusions.pdf \
    --annotation=!{gtf} \
    --cytobands=/opt2/arriba_v2.3.0/database/cytobands_hg19_hs37d5_GRCh37_v2.3.0.tsv \
    --proteinDomains=/opt2/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3
else
    touch !{dataset_id}.fusions.pdf
fi
    '''

}
