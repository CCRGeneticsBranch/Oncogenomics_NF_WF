include {DBinput_multiples} from '../modules/misc/DBinput'
include {Somatic_actionable} from '../modules/misc/DBinput'

workflow Annotation_somatic {


   somatic_actionable_sites = Channel.of(file(params.somatic_actionable_sites, checkIfExists:true))
   combined_gene_list = Channel.of(file(params.combined_gene_list, checkIfExists:true))
   group               = Channel.from("somatic")

take:
   dbinput_annot
   dbinput_somatic_snpeff
   dbinput_HC_snpeff
   dbinput_meta_tumor
   dbinput_meta_normal
   annotation_coding_rare

main:
   DBinput_multiples(
   dbinput_annot,
   dbinput_somatic_snpeff,
   dbinput_HC_snpeff,
   dbinput_meta_tumor,
   dbinput_meta_normal,
   group
   )

   Somatic_actionable(
   annotation_coding_rare,
   DBinput_multiples.out,
   somatic_actionable_sites,
   combined_gene_list

   )

emit:

somatic_actionable = Somatic_actionable.out
dbinput_somatic = DBinput_multiples.out

}
