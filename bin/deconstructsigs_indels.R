#!/usr/bin/env Rscript

library(BiocManager)
library(deconstructSigs)
library(BSgenome)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(colorspace)

args<-commandArgs(TRUE)
tsv_in <- args[1]
sample_name <- args[2]
cosmic_indel_rda <- args[3]
sigs.input = as.data.frame(t(read.table(tsv_in, sep="\t", header=TRUE, row.names=1)))
load(file= cosmic_indel_rda)
output.sigs = whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.indels.cosmic.v3.may2019, contexts.needed = TRUE)
sigs = data.frame(colnames(output.sigs$weights), t(output.sigs$weights))
sigs_file = paste(sample_name, "_cosmic_v3_data.csv", sep="")
write.csv(sigs, file=sigs_file)
bar_chart <- paste(sample_name, "_cosmic_v3.pdf", sep="")
pdf(bar_chart)
plotSignatures(output.sigs)
dev.off()
pie_chart <- paste(sample_name, "_cosmic_v3_pie.pdf", sep="")
num_sigs_shown <- sum(output.sigs$weights!=0)

all_sigs <- colnames(output.sigs$weights)

pie_colors <- sequential_hcl(17, h1 = 0, h2 = 360, c1 = 80, c2 = 100, l1 = 15, l2 = 90)
# pie_colors <- append(pie_colors, "#808080", after=length(pie_colors))
names(pie_colors) <- all_sigs
this_pie_pre <- pie_colors[output.sigs$weights!=0]
this_pie_pre_t <- t(as.data.frame(this_pie_pre))
pre_col <- colnames(this_pie_pre_t)
outweight_col <- colnames(output.sigs$weights)
this_pie <- this_pie_pre_t[,c(pre_col[pre_col %in% outweight_col])]
class(output.sigs)
pdf(pie_chart)

makePie(output.sigs, add.color = this_pie)
pie_labels = colnames(output.sigs$weights[c(output.sigs$weights!=0)])
pie_labels = append(pie_labels, "unknown", after=length(pie_labels))
this_pie <- append(this_pie, "#808080")
legend("bottomleft", legend=pie_labels, bty="n", fill=this_pie, cex = 0.8)
dev.off()
