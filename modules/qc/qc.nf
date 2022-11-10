process fastqc {
        tag { dataset_id }
        publishDir "$params.resultsdir/$dataset_id", mode: 'copy'

        input:
        tuple val(dataset_id),
        path(forward),
        path(reverse)

        output:
        tuple val("${dataset_id}"),
        path("fastqc_trim_${dataset_id}")


        script:
        """
        mkdir fastqc_trim_${dataset_id}
        fastqc -o fastqc_trim_${dataset_id} -q $forward $reverse
        """
}

process multiqc {
        tag { dataset_id }
        publishDir "$params.resultsdir/$dataset_id", mode: 'copy'
//        publishDir "s3://agc-424336837382-us-east-1/nfmvpout/$dataset_id", mode: 'copy'

        input:
	tuple val(dataset_id),
        path(qc)
        output:
        path "multiqc_report.html"

        script:
        """
	multiqc -m fastqc .

        """

}

process Genotyping {
        tag { dataset_id }

        input:
	tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome),
        path(Sites1000g4genotyping),
        path(vcf2genotype),
        path(vcf2loh)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}.star.samtools.vcf"),
        path("trim_${dataset_id}.star.gt"),
        path("trim_${dataset_id}.star.loh")

        script:
        """
	bcftools mpileup -O u -C50 $bam -f $genome -T $Sites1000g4genotyping   >trim_${dataset_id}.star.samtools.vcf

        perl $vcf2genotype trim_${dataset_id}.star.samtools.vcf >trim_${dataset_id}.star.gt

        perl $vcf2loh trim_${dataset_id}.star.samtools.vcf  > trim_${dataset_id}.star.loh
        """
}

