process Kraken {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc/kraken", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(meta), path(r1fq), path(r2fq),path(kraken_bacteria)

    output:
    tuple val(meta),
    path("${meta.lib}.kraken.taxa.txt"), emit: kraken_taxa
    tuple val(meta),
    path("${meta.lib}.krakenout"), emit : kraken_out
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}.kraken.taxa.txt"
    touch "${meta.lib}.krakenout"

    """


    script:
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    kraken --db ${kraken_bacteria} --fastq-input --gzip-compressed --threads ${task.cpus} --output ${prefix}.krakenout --preload --paired ${r1fq} ${r2fq}
    kraken-translate --mpa-format --db ${kraken_bacteria} ${prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > ${prefix}.krakentaxa

    mv ${prefix}.krakentaxa ${prefix}.kraken.taxa.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Kraken: \$(kraken --version|head -1|awk '{print \$3}')
    END_VERSIONS
    """
}

process Krona {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc/kraken", mode: 'copy'

    input:
    tuple val(meta), path(krakenout)

    output:
    tuple val(meta),
    path("${meta.lib}.krona.html")

    stub:
    """
    touch "${meta.lib}.krona.html"
    """

    script:
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    cut -f2,3 ${krakenout} |ktImportTaxonomy - -o ${prefix}.kronahtml
    mv ${prefix}.kronahtml ${prefix}.krona.html
    """
}





process Fastqc {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc/", mode: 'copy',pattern: "fastqc"

    input:
    tuple val(meta), path(samples)


    output:
    tuple val(meta), path("fastqc_${meta.lib}") , emit: fastqc_results
    path "versions.yml"             , emit: versions


    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    if [ ! -d fastqc_${meta.lib} ];then mkdir -p fastqc_${meta.lib};fi
    fastqc --extract ${samples.join(' ')} -t $task.cpus -o fastqc_${meta.lib}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Fastqc: \$(fastqc --version|awk '{print \$2}')
    END_VERSIONS
    """
}

process Fastq_screen {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc/fastq_screen/", mode: 'copy'

    input:
    tuple val(meta),path(trim),path(fastq_screen_config),path(fqs_db)

    stub:
    """
    touch "${meta.lib}_R1_screen.html"
    touch "${meta.lib}_R1_screen.html"
    """

    output:
    tuple val(meta),path("fastq_screen")


    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.lib}"

    """
    if [ ! -d fastq_screen ];then mkdir -p fastq_screen;fi
    ls ${fqs_db}
    fastq_screen --conf ${fastq_screen_config} --subset 1000000 --aligner bowtie2 --force ${trim[0]} ${trim[1]}  --outdir fastq_screen
    """
}


process Multiqc_old {
    tag "$meta.id"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}"

    input:
    path(qc_files)
    val(meta)


    output:
    tuple val(meta),path("multiqc_report.html")
    //path "versions.yml"

    script:
    """

    echo  "${qc_files.join('\n')}" > multiqc_input_files
    multiqc --file-list multiqc_input_files -f


    """

}

process Multiqc {
    tag "$meta.id"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}",pattern: "*html"

    input:
    tuple val(meta),path(qc_files)


    output:
    tuple val(meta), path("multiqc_report.html") , emit: multiqc_report
    path "versions.yml"             , emit: versions

    script:
    """

    echo  "${qc_files.join('\n')}" > multiqc_input_files

    multiqc --file-list multiqc_input_files -f

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    multiqc: \$(multiqc --version|sed 's/.*version //')
END_VERSIONS
    """

}

process Genotyping {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(Sites1000g4genotyping),
        path(genome),
        path(genome_fai),
        path(genome_dict)

    output:
    tuple val(meta),path("${meta.lib}.gt") , emit: gt
    tuple val(meta),path("${meta.lib}.loh"), emit: loh
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}.gt"
    touch "${meta.lib}.loh"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
     """
    bcftools mpileup -R ${Sites1000g4genotyping} -C50 -Ou -d 8000 -f ${genome} ${bam} | bcftools call --ploidy GRCh37 -m -Ov -o ${prefix}.samtools.vcf

    vcf2genotype.pl ${prefix}.samtools.vcf > ${prefix}.gt

    vcf2loh.pl ${prefix}.samtools.vcf  > ${prefix}.loh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version|grep -E '^bcftools'|sed 's/.*bcftools //')
    END_VERSIONS

    """
}

process Genotyping_Sample {

    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(gt_files)

    output:
    tuple val(meta),
        path("${meta.id}.genotyping.txt")

    stub:
    """
    touch "${meta.id}.genotyping.txt"
    """

    script:

     """
     mkdir -p GT RATIO
     cp ${gt_files.join(' ')} ./GT/
     echo Sample > RATIO/FirstColumn
     for file in GT/*
     do
                sample=`basename \$file .gt`
                echo \$sample
                echo \$sample >>RATIO/FirstColumn
                echo \$sample >>RATIO/\$sample.ratio
                for file2 in GT/*
                do
                        scoreGenotyes.pl \$file \$file2 >>RATIO/\$sample.ratio
                done
        done
     paste RATIO/FirstColumn RATIO/*.ratio >${meta.id}.genotyping.txt
     rm -rf GT RATIO
     sed -i 's/Sample_//g' ${meta.id}.genotyping.txt
     sed -i 's/.bwa//g' ${meta.id}.genotyping.txt
     sed -i 's/.star//g' ${meta.id}.genotyping.txt
     """


}


process CircosPlot_lib {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(loh)

    output:
    tuple val(meta),
        path("${meta.lib}.circos.png")

    stub:
    """
    touch "${meta.lib}.circos.png"
    """

    script:

     """
     circosLib.R  \$PWD/ ${meta.lib}.circos.png ${meta.lib}
     """
}

process CircosPlot {

    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(loh_files)

    output:
    tuple val(meta),
        path("${meta.id}.circos.png")

    stub:
    """
    touch "${meta.id}.circos.png"
    """

    script:

     """
     circos.R  \$PWD/ ${meta.id}.circos.png ${loh_files.join(' ')}
     """
}

process RNAseQC {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(rRNA_interval),
        path(transcript_gtf)

    output:
    tuple val(meta),
        path("rnaseqc")

    stub:
     """
     touch "report.html"
     """

    script:
     """
     java -jar \$RNASEQCJAR -r ${genome} -rRNA ${rRNA_interval} -o rnaseqc -s "${meta.lib}|${bam}|${meta.lib}" -t ${transcript_gtf}

     """

}

process Conpair_pile {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(conpair)

    output:
    tuple val(meta),
       path("${meta.lib}.conpair.mpileup")

    stub:
     """
       touch "${meta.lib}.conpair.mpileup"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     run_gatk_pileup_for_sample.py -B ${bam} -O ${prefix}.conpair.mpileup -R ${genome} -M ${conpair}
     """

}

process Exome_QC {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(Tbed),
        path(hsmetrics)

    output:
    tuple val(meta),path("${meta.lib}.consolidated_QC")

    stub:
     """
       touch "${meta.lib}.consolidated_QC"
     """

    script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     QC_stats_Final.py  ${bam} ${Tbed} . ${meta.id} ${prefix} "${meta.diagnosis}" ${hsmetrics} > ${prefix}.consolidated_QC
     """
}

process QC_summary_Patientlevel {

    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(qc_files)

    output:
    tuple val(meta),path("${meta.id}.consolidated_QC.txt")

    stub:
    """
      touch "${meta.id}.consolidated_QC.txt"
    """

    script:
    """
    cat ${qc_files.join(' ')} |grep '^#' |sort |uniq > header
    cat ${qc_files.join(' ')} |grep -v '^#' > sample_data
    cat header sample_data > ${meta.id}.consolidated_QC.txt

    """

}


process Conpair_concordance {

    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(tumor),
        path(normal),
        path(conpair)

    output:
    tuple val(meta),path("${meta.lib}.concod.txt")

    stub:
    """
      touch "${meta.lib}.concod.txt
    """

    script:
    """
    verify_concordance.py -T ${tumor} -N ${normal} -O ${meta.lib}.concod.txt  -C 20 -Q 30 -B 20 -M ${conpair}
    """

}

process Conpair_contamination {

    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(tumor),
        path(normal),
        path(conpair)

    output:
    tuple val(meta),path("${meta.lib}.conta.txt")

    stub:
    """
      touch "${meta.lib}.conta.txt
    """

    script:
    """
    estimate_tumor_normal_contamination.py -T ${tumor} -N ${normal} -O ${meta.lib}.conta.txt -Q 30 -M ${conpair}
    sed -i 's/Tumor sample contamination level: /${meta.lib}\t/g' ${meta.lib}.conta.txt
    sed -i 's/Normal sample contamination level: /${meta.normal_id}\t/g' ${meta.lib}.conta.txt
    """

}
