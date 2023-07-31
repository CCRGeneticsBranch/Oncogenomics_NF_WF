include {DBinput_multiples} from '../modules/misc/DBinput'
include {Germline_actionable} from '../modules/misc/DBinput'

workflow Annotation_germline {


   somatic_actionable_sites = Channel.of(file(params.somatic_actionable_sites, checkIfExists:true))
   combined_gene_list = Channel.of(file(params.combined_gene_list, checkIfExists:true))
   group               = Channel.from("germline")

take:
   dbinput_HC_annot
   dbinput_somatic_snpeff
   dbinput_HC_snpeff
   dbinput_meta_tumor
   dbinput_meta_normal
   annotation_coding_rare
   dbinput_somatic

main:
   DBinput_multiples(
   dbinput_HC_annot,
   dbinput_somatic_snpeff,
   dbinput_HC_snpeff,
   dbinput_meta_tumor,
   dbinput_meta_normal,
   group
   )

   Germline_actionable(
   annotation_coding_rare,
   DBinput_multiples.out,
   somatic_actionable_sites,
   combined_gene_list,
   dbinput_somatic

   )

emit:

germline = DBinput_multiples.out

}