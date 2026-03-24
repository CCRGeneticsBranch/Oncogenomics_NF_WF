#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(circlize))
library(stringr)

args<-commandArgs(TRUE)
DIR = str_trim(args[1])
FILE=str_trim(args[2])
GENOME=str_trim(args[4])

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


png(FILE, width=10, height=10, units="in", res=200, type=c("cairo"))

# Initialize circos with specified genome (hg19 or hg38)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species=GENOME)

# Create a track for each sample
for (i in 1:length(files)){
        LOH.data <- read.table(paste(DIR, files[i], sep = ""), sep="\t", quote="", header=TRUE)

        # Create empty track with fixed scale 0-1
        circos.track(ylim=c(0, 1))

        # Plot all points for this sample across all chromosomes (simpler & faster)
        circos.trackPoints(LOH.data[,1], LOH.data[,2], LOH.data[,3],
                          col = cols[i], pch = 16, cex = 0.3)
}

# Clear circos
circos.clear()

# Add legend
legend("topleft", legend=labs, col=cols, pch=19, cex=0.80, bty="n")

dev.off()
