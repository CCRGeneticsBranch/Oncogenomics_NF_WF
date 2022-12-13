
workflow Mixcr_VCJtools{
    take:
        trimmed_fq_w_mixcr_license
    main:
        Mixcr(trimmed_fq_w_mixcr_license)
        VDJtools(Mixcr.out)
    emit:
        VDJtools.out
}

process Mixcr{
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/mixcr", mode: "copy"

    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2),
        path(license)
    
    output:
    tuple val("${dataset_id}"),
        path("mixcr")


    shell:
    '''
set -exo pipefail
if [ -d /lscratch/${SLURM_JOB_ID} ];then
    TMPDIR="/lscratch/${SLURM_JOB_ID}/!{dataset_id}_mixcr"
elif [ -d /dev/shm ];then
    TMPDIR="/dev/shm/!{dataset_id}_mixcr"
else 
    TMPDIR="/tmp/!{dataset_id}_mixcr"
fi
if [ -d ${TMPDIR} ];then rm -rf ${TMPDIR}; fi

echo ""
echo "Mixcr LICENSE"
license_path=`pwd`
echo $license_path
export MI_LICENSE_FILE="$license_path/mi.license"
echo $MI_LICENSE_FILE
cat $MI_LICENSE_FILE
echo ""


if [ -d mixcr ];then rm -rf mixcr;fi
mkdir -p mixcr && cd mixcr

# using 4.1.0 mixr version

mixcr analyze rnaseq-tcr-full-length --species hsa --rna ../!{r1} ../!{r2} !{dataset_id}

if [ -f "!{dataset_id}.clones.tsv" ];then 
    mv "!{dataset_id}.clones.tsv" "!{dataset_id}.clones.RNA.txt"
else
    touch "!{dataset_id}.clones.RNA.txt"
fi

mixcr exportClones -c TRA !{dataset_id}.contigs.clns !{dataset_id}.clones.TRA.txt
mixcr exportClones -c TRB !{dataset_id}.contigs.clns !{dataset_id}.clones.TRB.txt
mixcr exportClones -c IGH !{dataset_id}.contigs.clns !{dataset_id}.clones.IGH.txt
mixcr exportClones -c IGK !{dataset_id}.contigs.clns !{dataset_id}.clones.IGK.txt
mixcr exportClones -c IGL !{dataset_id}.contigs.clns !{dataset_id}.clones.IGL.txt

for f in TRA TRB IGH IGK IGL;do
    if [ ! -f "!{dataset_id}.clones.${f}.txt" ];then
        touch "!{dataset_id}.clones.${f}.txt"
    fi
done

#Older mixcr version

#mixcr align -t !{task.cpus} -p rnaseq-cdr3 -s hsa -OallowPartialAlignments=true --report !{dataset_id}.alignment.log ../!{r1} ../!{r2} !{dataset_id}.aln.vdjca

#mixcr assemblePartial !{dataset_id}.aln.vdjca alignments_rescued_1.vdjca
#mixcr assemblePartial alignments_rescued_1.vdjca alignments_rescued_2.vdjca

#mixcr extendAlignments alignments_rescued_2.vdjca alignments_rescued_2_extended.vdjca

#mixcr assemble -t !{task.cpus} alignments_rescued_2_extended.vdjca !{dataset_id}.clones.clns

#mixcr exportClones !{dataset_id}.clones.clns !{dataset_id}.clones.RNA.txt
#mixcr exportClones -c TRA !{dataset_id}.clones.clns !{dataset_id}.clones.TRA.txt
#mixcr exportClones -c TRB !{dataset_id}.clones.clns !{dataset_id}.clones.TRB.txt
#mixcr exportClones -c IGH !{dataset_id}.clones.clns !{dataset_id}.clones.IGH.txt
#mixcr exportClones -c IGK !{dataset_id}.clones.clns !{dataset_id}.clones.IGK.txt
#mixcr exportClones -c IGL !{dataset_id}.clones.clns !{dataset_id}.clones.IGL.txt

cd ..
'''
}


process VDJtools{
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/mixcr", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(mixcrDir)
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}.summarystats.RNA.txt")


    shell:
    '''
set -exo pipefail
if [ -d /lscratch/${SLURM_JOB_ID} ];then
    TMPDIR="/lscratch/${SLURM_JOB_ID}/!{dataset_id}_VDJtools"
elif [ -d /dev/shm ];then
    TMPDIR="/dev/shm/!{dataset_id}_VDJtools"
else 
    TMPDIR="/tmp/!{dataset_id}_VDJtools"
fi
if [ -d ${TMPDIR} ];then rm -rf ${TMPDIR}; fi

cd mixcr

java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} Convert -S mixcr !{dataset_id}.clones.RNA.txt convert
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} Convert -S mixcr !{dataset_id}.clones.TRA.txt convert
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} Convert -S mixcr !{dataset_id}.clones.TRB.txt convert
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} Convert -S mixcr !{dataset_id}.clones.IGH.txt convert
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} Convert -S mixcr !{dataset_id}.clones.IGK.txt convert
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} Convert -S mixcr !{dataset_id}.clones.IGL.txt convert

java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} CalcBasicStats convert.!{dataset_id}.clones.RNA.txt !{dataset_id}.RNA
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} CalcBasicStats convert.!{dataset_id}.clones.TRA.txt !{dataset_id}.TRA
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} CalcBasicStats convert.!{dataset_id}.clones.TRB.txt !{dataset_id}.TRB
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} CalcBasicStats convert.!{dataset_id}.clones.IGH.txt !{dataset_id}.IGH
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} CalcBasicStats convert.!{dataset_id}.clones.IGK.txt !{dataset_id}.IGK
java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} CalcBasicStats convert.!{dataset_id}.clones.IGL.txt !{dataset_id}.IGL

# getting this error hence turning ploting off for now
# Command output:
#   ========================================
#   circlize version 0.4.15
#   CRAN page: https://cran.r-project.org/package=circlize
#   Github page: https://github.com/jokergoo/circlize
#   Documentation: https://jokergoo.github.io/circlize_book/book/
#   
#   If you use it in published research, please cite:
#   Gu, Z. circlize implements and enhances circular visualization
#     in R. Bioinformatics 2014.
#   
#   This message can be suppressed by:
#     suppressPackageStartupMessages(library(circlize))
#   ========================================
#   
#   Loading required package: RColorBrewer
#   Error in apply(temp[1, 2:m], 2, as.character) : 
#     dim(X) must have a positive length
#   Execution halted
# java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} PlotFancyVJUsage convert.!{dataset_id}.clones.RNA.txt !{dataset_id}.RNA
# java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} PlotFancyVJUsage convert.!{dataset_id}.clones.TRA.txt !{dataset_id}.TRA
# java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} PlotFancyVJUsage convert.!{dataset_id}.clones.TRB.txt !{dataset_id}.TRB
# java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} PlotFancyVJUsage convert.!{dataset_id}.clones.IGH.txt !{dataset_id}.IGH
# java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} PlotFancyVJUsage convert.!{dataset_id}.clones.IGK.txt !{dataset_id}.IGK
# java -Djava.io.tmpdir=${TMPDIR} -jar ${VDJTOOLS_JAR} PlotFancyVJUsage convert.!{dataset_id}.clones.IGL.txt !{dataset_id}.IGL

head -n1 !{dataset_id}.RNA.basicstats.txt > !{dataset_id}.summarystats.RNA.txt
cat !{dataset_id}.TRA.basicstats.txt !{dataset_id}.TRB.basicstats.txt !{dataset_id}.IGH.basicstats.txt !{dataset_id}.IGK.basicstats.txt !{dataset_id}.IGL.basicstats.txt | grep -v sample_id >> !{dataset_id}.summarystats.RNA.txt

cd ..
cp mixcr/!{dataset_id}.summarystats.RNA.txt .
ls -alrth

'''

}
