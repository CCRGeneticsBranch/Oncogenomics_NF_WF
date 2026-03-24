#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(circlize))
library(stringr)

args<-commandArgs(TRUE)
DIR = str_trim(args[1])
FILE = str_trim(args[2])
GENOME = str_trim(args[3])
SAM = str_trim(args[4])


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

# Initialize circos with specified genome
circos.par("start.degree" = 90, "track.height" = 0.1, "gap.degree" = 2)
circos.initializeWithIdeogram(species=GENOME)


z = length(files)

# Determine track height and point size based on number of samples
if(z >= 4){
    # Many samples: use available space efficiently
    track_height = 0.08
    point_cex = 0.15
} else if(z == 1) {
    # Single sample: larger track for better visibility
    track_height = 0.15
    point_cex = 0.3
} else {
    # 2-3 samples: medium sized tracks
    track_height = 0.10
    point_cex = 0.25
}

# Create a track for each sample
for (i in 1:length(files)){
    LOH.data <- read.table(paste(DIR, files[i], sep = ""), sep="\t", quote="", header=TRUE)

    # Create empty track with fixed scale 0-1 and dynamic height
    circos.track(ylim=c(0, 1), track.height=track_height, bg.border=NA)

    # Plot all points for this sample across all chromosomes
    # Use dynamically sized points based on sample count
    circos.trackPoints(LOH.data[,1], LOH.data[,2], LOH.data[,3],
                      col = cols[i], pch = 16, cex = point_cex)
}

# Clear circos
circos.clear()

# Dynamic legend placement based on number of samples
x = round(length(files)/2)
y = x + 1

if(z > 8){
    # Split legend for many samples
    legend("topright", legend=labs[1:x], col=cols[1:x], pch=19, cex=0.80, bty="n")
    legend("topleft", legend=labs[y:z], col=cols[y:z], pch=19, cex=0.80, bty="n")
} else {
    # Single legend for few samples
    legend("topleft", legend=labs[1:z], col=cols[1:z], pch=19, cex=0.80, bty="n")
}

dev.off()
