process Manta {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(Tbam),path(Tindex),path(Tbed)
     tuple val(meta2),path(Nbam),path(Nindex)
     path genome
     path genome_fai
     path genome_dict
     
     output:
     tuple val(meta),path("manta")

    stub:
    """   
        touch "candidateSmallIndels.vcf.gz"

    """

    script:

    """
    configManta.py --normalBam ${Nbam} --tumorBam ${Tbam} --referenceFasta ${genome} --runDir ./manta
    ./manta/runWorkflow.py
    """
}
