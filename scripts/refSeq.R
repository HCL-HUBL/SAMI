#!/usr/bin/env Rscript

### Prepare introns and exons tracks from a refGene file from UCSC
### Author : <sylvain.mareschal@chu-lyon.fr>

# Collect CLI arguments
args <- commandArgs(TRUE)
if(length(args) != 4L) stop("USAGE : ./refSeq.R refGene.txt ORGANISM ASSEMBLY CHROMOSOMES")
file <- args[1]
organism <- args[2]
assembly <- args[3]
chromosomes <- strsplit(args[4], split=",")[[1]]

# Check CLI arguments
if(!file.exists(file)) stop("Provided refGene file must exist")

### file <- "store/gencode.v35.primary_assembly.annotation.gtf"
### organism <- "Human"
### assembly <- "GRCh38"
### chromosomes <- c(1:22, "X", "Y")

library(Rgb)

# UCSC refGene
tab <- read.table(
	file, sep="\t", stringsAsFactors = FALSE,
	col.names = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"),
)

# Transcripts annotated but of interest (EGFR)
tab <- tab[ tab$name != "NM_001346941" ,]

# Collect exon junctions
# Mostly 'chrom:left-right', very few exceptions on unusual chromosomes and left = right + 1
starts <- strsplit(tab$exonStarts, split=",", fixed=TRUE)
ends <- strsplit(tab$exonEnds, split=",", fixed=TRUE)
junctions <- list()
for(i in 1:nrow(tab)) junctions[[i]] <- sprintf("%s:%i-%s", tab$chrom[i], as.integer(head(ends[[i]], -1)) + 1L, tail(starts[[i]], -1))
jun <- sub("^chr", "", unique(unlist(junctions)))
saveRDS(jun, file=sprintf("introns.refSeq.%s.rds", assembly))

# Collect exons
trk <- list()
for(i in 1:nrow(tab)) {
	if(i %% 1000L == 0L) message(i, "/", nrow(tab))
	trk[[i]] <- data.frame(
		name       = "",
		chrom      = sub("^chr", "", tab$chrom[i]),
		strand     = tab$strand[i],
		start      = as.integer(starts[[i]]),
		end        = as.integer(ends[[i]]),
		transcript = sprintf("%s (%i)", tab$name[i], i),
		symbol     = tab$name2[i],
		stringsAsFactors = FALSE
	)
}
trk <- do.call(rbind, trk)

# Transcript track
trk <- trk[ trk$chrom %in% chromosomes ,]
trk$chrom <- factor(trk$chrom, levels=chromosomes)
trk$strand <- factor(trk$strand, levels=c("-","+"))
trk <- track.exons(
	trk,
	.name = "RefSeq",
	.organism = organism,
	.assembly = assembly
)
trk$buildGroupPosition("transcript")
trk$buildGroupSize("transcript")
trk$setParam("maxElements", 1000)

# Export
saveRDT(trk, file=sprintf("exons.refSeq.%s.rdt", assembly))
