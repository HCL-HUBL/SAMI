#!/usr/bin/env Rscript

### Converts a GENCODE GTF to the refFlat format required by GATK CollectRnaSeqMetrics
### Author : <sylvain.mareschal@lysarc.org>

# CLI arguments
args <- commandArgs(TRUE)
if(length(args) != 2L) stop("USAGE: ./gtfToReflat.R INPUT.gtf OUTPUT.refFlat")
inputFile <- args[1]
outputFile <- args[2]

# Check arguments
if(!file.exists(inputFile)) stop("'INPUT.gtf' must exist")
if(file.exists(outputFile)) warning("'OUTPUT.refFlat' will be replaced")

# GTF transcriptome definition to parse
con <- file(inputFile, open="rt")
fileSize <- as.double(file.info(inputFile)["size"])

# Storage
trs <- trs.chrom <- trs.strand <- trs.start <- trs.end <- trs.gene <- NULL
cds <- cds.start <- cds.end <- NULL
exons <- exons.start <- exons.end <- NULL

# Parse by chunk
chunkSize <- 30000L
repeat{
	# Parse GTF/GFF3 file
	gtf <- read.table(
		file = con, sep = "\t", header = FALSE, quote = NULL, nrow=chunkSize,
		col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
		colClasses = c("character", "character", "character", "integer", "integer", "character", "character", "character", "character"),
	)
	if(nrow(gtf) == 0L) break

	# Get transcript definition
	filter <- gtf$feature == "transcript"
	trs <- c(trs, sub("^.+; transcript_id \"([^\"]+)\";.+$", "\\1", gtf$attribute[filter]))
	trs.chrom <- c(trs.chrom, gtf$seqname[filter])
	trs.strand <- c(trs.strand, gtf$strand[filter])
	trs.start <- c(trs.start, gtf$start[filter])
	trs.end <- c(trs.end, gtf$end[filter])
	trs.gene <- c(trs.gene, sub("^.+; gene_name \"([^\"]+)\";.+$", "\\1", gtf$attribute[filter]))
	
	# Get CDS boundaries
	filter <- gtf$feature == "CDS"
	cds <- c(cds, sub("^.+; transcript_id \"([^\"]+)\";.+$", "\\1", gtf$attribute[filter]))
	cds.start <- c(cds.start, gtf$start[filter])
	cds.end <- c(cds.end, gtf$end[filter])

	# Get exon boundaries
	filter <- gtf$feature == "exon"
	exons <- c(exons, sub("^.+; transcript_id \"([^\"]+)\";.+$", "\\1", gtf$attribute[filter]))
	exons.start <- c(exons.start, gtf$start[filter])
	exons.end <- c(exons.end, gtf$end[filter])
	
	# Progression (requires file size)
	message(length(trs), " transcripts processed (", round(100 * seek(con) / fileSize, 1), "%)")
}

# Close file
close(con)

# Get CDS starts and ends
cdsStarts <- tapply(X=cds.start, INDEX=cds, FUN=min)
cdsEnds <- tapply(X=cds.end, INDEX=cds, FUN=max)

# Group exon coordinates
exonCounts <- tapply(X=exons.start, INDEX=exons, FUN=length)
exonStarts <- tapply(X=exons.start, INDEX=exons, FUN=function(x) paste(sort(x), collapse=","))
exonEnds <- tapply(X=exons.end, INDEX=exons, FUN=function(x) paste(sort(x), collapse=","))

# Reshape as refFlat
# <http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?hgta_doSchemaDb=hg19&hgta_doSchemaTable=refFlat>
out <- data.frame(
	geneName = trs.gene,
	name = trs,
	chrom = trs.chrom,
	strand = trs.strand,
	txStart = trs.start,
	txEnd = trs.end,
	stringsAsFactors = FALSE
)

# CDS
i <- match(trs, names(cdsStarts))
out$cdsStart <- as.integer(cdsStarts[i])
out$cdsEnd <- as.integer(cdsEnds[i])

# Exons
i <- match(trs, names(exonCounts))
out$exonCounts <- as.integer(exonCounts[i])
out$exonStarts <- as.character(exonStarts[i])
out$exonsEnds <- as.character(exonEnds[i])

# Fill NAs for non-coding transcripts as in UCSC files
# <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz>
out[ is.na(out$cdsStart) , "cdsStart" ] <- out[ is.na(out$cdsStart) , "txEnd" ]
out[ is.na(out$cdsEnd) , "cdsEnd" ] <- out[ is.na(out$cdsEnd) , "txEnd" ]

# Add extra comma to exon lists as in UCSC files
# <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz>
out$exonStarts <- sprintf("%s,", out$exonStarts)
out$exonsEnds <- sprintf("%s,", out$exonsEnds)

# Checks
if(!identical(names(cdsStarts), names(cdsEnds)))     stop("CDS index inconsistency")
if(!identical(names(exonCounts), names(exonStarts))) stop("exon index inconsistency (start)")
if(!identical(names(exonCounts), names(exonEnds)))   stop("Exon index inconsistency (end)")
if(any(is.na(out$exonCounts) | out$exonCounts < 1L)) stop("Transcript without exon")
if(any(!cds %in% trs))                               stop("Unknown CDS")
if(any(!exons %in% trs))                             stop("Unknown Exon")

# Export as TSV to stdout
write.table(out, file=outputFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, na="")

