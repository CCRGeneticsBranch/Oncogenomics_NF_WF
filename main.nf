params.reads = "s3://myncitestbucket/inputs/ggal_gut_{1,2}.fastq"
reads_ch = Channel.fromFilePairs(params.reads)

process cutadapt {
	container 'nciccbr/ncigb_cutadapt_v1.18:latest'
	publishDir 's3://agc-913060503860-us-west-2/nfmvpout', mode: 'copy'
	
	input:
	tuple val(sample_id), path(reads) from reads_ch

	output:
	path "trim*" into trim_ch1, trim_ch2

	
	"""
	cutadapt  -o trim_${sample_id}_R1.fastq -p trim_${sample_id}_R2.fastq ${reads[0]} ${reads[1]}
	"""

}
