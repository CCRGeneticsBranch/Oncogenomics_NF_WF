#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(OmicCircos))
library(stringr)
library(RColorBrewer)

args<-commandArgs(TRUE)
DIR = str_trim(args[1])
FILE=str_trim(args[2])
SAM=str_trim(args[3])

files <- list.files(path = DIR, pattern=".loh$")

labs <- paste("", gsub("Sample_|\\.bwa|\\.star|\\.loh", "", files, perl=TRUE), sep="")


cols <-c('#26294a','#01545a','#bd544f','#017351',
	'#03c383','#b8bd4f','#aad962','#fbbf45',
	'#bd8b4f','#ef6a32','#ed0346','#d76e60',
	'#a12a5e','#710162','#26294a','#01545a',
	'#bd544f','#017351','#03c383','#b8bd4f',
	'#aad962','#fbbf45','#bd8b4f','#ef6a32',
	'#ed0346','#d76e60','#a12a5e','#710162'
       )

options(stringsAsFactors = FALSE)
set.seed(1234)
png(FILE ,width = 2000, height = 2000, res=200, points=12, type=c("cairo"))
par(mar=c(2, 2, 4, 2))
plot(c(1,900), c(1,900), type="n", axes=FALSE, xlab="", ylab="", main="")

circos(R=400, cir="hg19", type="chr", mapping=UCSC.hg19.chr,print.chr.lab=TRUE, W=10, lwd=5, cex=1.5)

r=300

for (i in 1:length(files)){
        LOH.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=T)
        circos(cir="hg19", R=r, W=100, type="s", mapping=LOH.data, col.v=3, col=cols[i], B=FALSE, cex=0.0003, lwd=1)
	r=r-100
}



legend("topleft", legend=labs, col=cols,pch=19,cex= 0.80,bty="n")
legendstyle =list("x" = 0.5, "y" = -100)
layoutstyle = list(legend=legendstyle)
dev.off()
