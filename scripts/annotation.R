#!/usr/bin/env Rscript

### Prepare introns and exons tracks from a GTF file
### Author : <sylvain.mareschal@chu-lyon.fr>

# Collect CLI arguments
args <- commandArgs(TRUE)
if(length(args) != 4L) stop("USAGE : ./annotation.R ANNOTATION.gtf ORGANISM ASSEMBLY CHROMOSOMES")
file <- args[1]
organism <- args[2]
assembly <- args[3]
chromosomes <- strsplit(args[4], split=",")[[1]]

# Check CLI arguments
if(!file.exists(file)) stop("ANNOTATION.gtf must exist")

library(Rgb)

# Parse exons
gtf <- read.gtf(pipe(sprintf("awk '$3 == \"exon\" { print }' \"%s\"", file), "rt"))

# Black list
exclude <- c(
	grep("^NM_001346941", gtf$transcript_id) ### EGFR pathological transcript
)
if(length(exclude) > 0L) gtf <- gtf[ -exclude ,]

# Possible symbol columns, ordered by priority
columns <- c("gene_name", "gene", "gene_id")
column <- columns[ columns %in% colnames(gtf) ][1]
if(is.na(column)) stop("No suitable column found for symbol")

# Group and sort by transcript
chrom <- tapply(X=gtf$seqname, INDEX=gtf$transcript_id, FUN=unique)
strand <- tapply(X=as.character(gtf$strand), INDEX=gtf$transcript_id, FUN=unique)
symbol <- tapply(X=gtf[[column]], INDEX=gtf$transcript_id, FUN=unique)
starts <- tapply(X=gtf$start, INDEX=gtf$transcript_id, FUN=sort)
ends <- tapply(X=gtf$end, INDEX=gtf$transcript_id, FUN=sort)

# Collect exon junctions as a 'chrom:left-right' vector
junctions <- list()
for(i in 1:length(starts)) junctions[[i]] <- sprintf("%s:%i-%s", chrom[i], as.integer(head(ends[[i]], -1)) + 1L, tail(starts[[i]], -1))
jun <- sub("^chr", "", unique(unlist(junctions)))
saveRDS(jun, file=sprintf("introns.%s.rds", assembly))

# Collect exons
n.exons <- sapply(starts, length)
tab <- data.frame(
	name       = "",
	chrom      = rep(sub("^chr", "", chrom), n.exons),
	strand     = rep(strand, n.exons),
	start      = as.integer(unlist(starts)),
	end        = as.integer(unlist(ends)),
	transcript = rep(names(chrom), n.exons),
	symbol     = rep(symbol, n.exons),
	stringsAsFactors = FALSE
)

# Transcript track
tab <- tab[ tab$chrom %in% chromosomes ,]
tab$chrom <- factor(tab$chrom, levels=chromosomes)
tab$strand <- factor(tab$strand, levels=c("-","+"))
trk <- track.exons(
	tab,
	.name = basename(file),
	.organism = organism,
	.assembly = assembly
)
trk$buildGroupPosition("transcript")
trk$buildGroupSize("transcript")
trk$setParam("maxElements", 1000)

# Export
saveRDT(trk, file=sprintf("exons.%s.rdt", assembly))

