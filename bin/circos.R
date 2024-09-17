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
par(mar=c(0, 1, 0, 1))
plot(c(1,900), c(1,900), type="n", axes=FALSE, xlab="", ylab="", main="")

#seg.len= (r-100)/length(files)
circos(R=380, cir="hg19", type="chr", mapping=UCSC.hg19.chr,print.chr.lab=TRUE, W=10, lwd=5, cex=1.5)

r.start=380
seg.len= (r.start-100)/length(files)
r = r.start - seg.len;

z=length(files)
if(z >= 4){
         for (i in 1:length(files)){
               LOH.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=T)
               circos(cir="hg19", R=r, W=seg.len-0.5, type="s", mapping=LOH.data, col.v=3, col=cols[i], B=FALSE, cex=0.0003, lwd=1)
               r=r-seg.len 
         }
} else if(z == 1) {
         seg.len= (r.start-325)/length(files)
         r = r.start - seg.len;
         for (i in 1:length(files)){
               LOH.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=T)
               circos(cir="hg19", R=r, W=seg.len-0.5, type="s", mapping=LOH.data, col.v=3, col=cols[i], B=FALSE, cex=0.0003, lwd=1)
               r=r-seg.len
         }
} else {
         seg.len= (r.start-275)/length(files)
         r = r.start - seg.len;
         for (i in 1:length(files)){
               LOH.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=T)
               circos(cir="hg19", R=r, W=seg.len-0.5, type="s", mapping=LOH.data, col.v=3, col=cols[i], B=FALSE, cex=0.0003, lwd=1)
               r=r-seg.len
         }
}


x=length(files)/2
x=round(x)
y= x+1
z=length(files)


if(z > 8){
        legend("topright", legend=labs[1:x], col=cols[1:x],pch=19,cex= 0.80,bty="n")
        legend("topleft", legend=labs[y:z], col=cols[y:z],pch=19,cex=0.80,bty="n")
} else {
        legend("topleft", legend=labs[1:z], col=cols[1:z],pch=19,cex= 0.80,bty="n")
}

###legend("bottom", legend=labs, col=cols, lty=1, lwd=4, cex=.75, xpd = TRUE,inset=c(0,-0.25))

legendstyle =list("x" = 0.5, "y" = -100)
layoutstyle = list(legend=legendstyle)
dev.off()
