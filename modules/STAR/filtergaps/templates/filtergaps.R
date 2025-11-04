#!/usr/bin/env Rscript --vanilla

# Collect Nextflow arguments
junctionFile <- "!{junctions}"
depthFile <- "!{depth}"
threshold <- as.double("!{threshold}")



message("Parsing junctions...")

junctions <- read.table(
	junctionFile, sep="\t", quote=NULL, comment.char="",
	col.names = c("chrom", "start", "end", "strand", "motif", "annotated", "reads.uni", "reads.multi", "overhang"),
	colClasses = c("character", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer")
)

message("Preparing positions...")

positions <- rbind(
	data.frame(chrom=junctions$chrom, pos=junctions$start-1L),
	data.frame(chrom=junctions$chrom, pos=junctions$end+1L)
)
positions <- positions[ order(positions$chrom, positions$pos) ,]
positions <- unique(positions)
write.table(positions, file="positions.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

message("Querying depth...")

system(sprintf("tabix -R \"positions.tsv\" \"%s\" > positions.out", depthFile), intern=TRUE)

message("Reshapping depth results...")

tmp <- read.table(
	"positions.out", sep="\t", quote=NULL, comment.char="",
	col.names = c("chrom", "pos", "depth"),
	colClasses = c("character", "integer", "integer")
)
depth <- tmp$depth
names(depth) <- paste(tmp$chrom, tmp$pos, sep=":")

message("Merging...")

junctions$start.depth <- depth[ paste(junctions$chrom, junctions$start-1L, sep=":") ]
junctions$end.depth <- depth[ paste(junctions$chrom, junctions$end+1L, sep=":") ]
junctions$start.depth[ is.na(junctions$start.depth) ] <- 0L
junctions$end.depth[ is.na(junctions$end.depth) ] <- 0L

message("Computing score...")

junctions$score <- (junctions$reads.uni + junctions$reads.multi) / pmax(junctions$start.depth, junctions$end.depth)

message("Filtering...")

out <- junctions[ junctions$score >= threshold , 1:9 ]

message("Exporting...")

dir.create("out")
write.table(out, file=sprintf("out/%s", junctionFile), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(junctions, file="scored-junctions.tsv", sep="\t", row.names=FALSE)
