#!/usr/bin/env Rscript --vanilla

# Collect Nextflow arguments
junctionFile <- "!{junctions}"
depthFile <- "!{depth}"
threshold <- "!{threshold}"



library(Rgb)

# Parse file
tab <- read.table(
	junctionFile, sep="\t", quote=NULL, comment.char="",
	col.names = c("chrom", "start", "end", "strand", "motif", "annotated", "reads.uni", "reads.multi", "overhang"),
	colClasses = c("character", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer"),
)

# Storage
score <- double(nrow(tab))

for(i in 1:nrow(tab)) {
    # Shortcuts
	chrom <- tab$chrom[i]
    start <- tab$start[i]
    end <- tab$end[i]

    # Depth at starting site
    depth.start <- as.integer(system(sprintf("tabix \"%s\" %s:%i-%i | cut -f3", depthFile, chrom, start-1L, start-1L), intern=TRUE))
    if(length(depth.start) == 0L) depth.start <- 0L

    # Depth at ending site
    depth.end <- as.integer(system(sprintf("tabix \"%s\" %s:%i-%i | cut -f3", depthFile, chrom, end+1L, end+1L), intern=TRUE))
    if(length(depth.end) == 0L) depth.end <- 0L

    # Reads supporting the junction
    reads <- tab$reads.uni[i] + tab$reads.multi[i]
	
	# Normalized support
    score[i] <- reads / max(depth.start, depth.start)
}

# Filter out junctions with low normalized support
tab <- tab[ score >= threshold ,]

# Export
dir.create("out")
write.table(tab, file=sprintf("out/%s", junctionFile), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
