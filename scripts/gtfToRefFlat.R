#!/usr/bin/env Rscript

### Converts a GENCODE GTF to the refFlat format required by GATK CollectRnaSeqMetrics
### <http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?hgta_doSchemaDb=hg19&hgta_doSchemaTable=refFlat>
### GTF have 1-based start & end, UCSC uses 0-based start : <http://genome.ucsc.edu/FAQ/FAQtracks#tracks1>

# CLI arguments
args <- commandArgs(TRUE)
if(length(args) != 2L) stop("USAGE: ./gtfToReflat.R INPUT.gtf OUTPUT.refFlat")
inputFile <- args[1]
outputFile <- args[2]

# Check arguments
if(!file.exists(inputFile)) stop("'INPUT.gtf' must exist")
if(file.exists(outputFile)) warning("'OUTPUT.refFlat' will be replaced")

library(Rgb)

# Parse GTF
CDS <- read.gtf(pipe(sprintf("awk '$3 == \"CDS\" { print }' \"%s\"", inputFile), "rt"))
exons <- read.gtf(pipe(sprintf("awk '$3 == \"exon\" { print }' \"%s\"", inputFile), "rt"))
transcript <- read.gtf(pipe(sprintf("awk '$3 == \"transcript\" { print }' \"%s\"", inputFile), "rt"))

# Possible symbol columns, ordered by priority
columns <- c("gene_name", "gene", "gene_id")
column <- columns[ columns %in% colnames(transcript) ][1]
if(is.na(column)) stop("No suitable column found for symbol")

# Transcript-based
out <- data.frame(
	"geneName" = transcript[[ column ]],
	"name"     = transcript$transcript_id,
	"chrom"    = transcript$seqname,
	"strand"   = transcript$strand,
	"txStart"  = transcript$start - 1L,
	"txEnd"    = transcript$end,
	stringsAsFactors = FALSE
)

# Add CDS boundaries (one row per exon)
cds.start <- tapply(X=CDS$start - 1L, INDEX=CDS$transcript_id, FUN=min)
cds.end <- tapply(X=CDS$end, INDEX=CDS$transcript_id, FUN=max)
out$"cdsStart" <- cds.start[ out$name ]
out$"cdsEnd" <- cds.end[ out$name ]

# Add exons
starts <- tapply(X=exons$start - 1L, INDEX=exons$transcript_id, FUN=sort)
ends <- tapply(X=exons$end, INDEX=exons$transcript_id, FUN=sort)
out$"exonCount" <- sapply(starts, length)[ out$name ]
out$"exonStarts" <- sapply(starts, paste, collapse=",")[ out$name ]
out$"exonEnds" <- sapply(ends, paste, collapse=",")[ out$name ]

# Fill NAs for non-coding transcripts as in UCSC files
# <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz>
out[ is.na(out$cdsStart) , "cdsStart" ] <- out[ is.na(out$cdsStart) , "txEnd" ]
out[ is.na(out$cdsEnd) , "cdsEnd" ] <- out[ is.na(out$cdsEnd) , "txEnd" ]

# Add extra comma to exon lists as in UCSC files
# <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz>
out$exonStarts <- sprintf("%s,", out$exonStarts)
out$exonEnds <- sprintf("%s,", out$exonEnds)

# Check
if(any(is.na(out$exonCounts) | out$exonCounts < 1L)) stop("Transcript without exon")

# Export as TSV to stdout
write.table(out, file=outputFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, na="")

