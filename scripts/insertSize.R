#!/usr/bin/env Rscript --vanilla

### Estimates the distribution of insert sizes from the first read pairs in a BAM file (use a **transcriptome-aligned** BAM file for RNA-seq !)
### Prints to stdout the YAML formatted results, for handling by MultiQC (automatically if stored in a file with name "*_mqc.yaml")
### Requires SAMtools to be available from the PATH
### 
### EXAMPLE : ./insertSize.R sampleA ../STAR/sampleA/Aligned.toTranscriptome.out.bam > sampleA_mqc.yaml
### EXAMPLE : ./insertSize.R sampleA ../STAR/sampleA/isize.txt > sampleA_mqc.yaml
### 
### Author : <sylvain.mareschal@lysarc.org>

nReads <- 1e6
nPoints <- 100L



# CLI arguments
args <- commandArgs(TRUE)
if(!length(args) %in% 2:3) stop("USAGE: ./insertSize.R SAMPLE_ID INPUT_FILE [ SAMTOOLS ]")
sample <- args[1]
inputFile <- args[2]
samtools <- args[3]

# Check arguments
if(samtools == "" || is.na(samtools)) samtools <- "samtools"
if(!file.exists(inputFile)) stop("INPUT_FILE must exist")



if(grepl("\\.bam", inputFile, ignore.case=TRUE)) {
	# BAM input, get ISIZE field for nth first alignments (R2 & proper-paired only)
	isize <- system(sprintf("\"%s\" view -f 0x2 -f 0x80 \"%s\" | cut -f9 | head -%i", samtools, inputFile, nReads), intern=TRUE)
	isize <- abs(as.integer(isize))
} else {
	# Assume text input corresponding to the command above
	isize <- abs(scan(inputFile, what=0L, sep="\n", quiet=TRUE))
}

if(length(isize) == 0L) {
	# No proper pair
	x <- 0:100
	y <- rep(0, length(x))
} else {
	# Histogram
	binSize <- 10L
	breaks <- seq(
		from = binSize * (floor(min(isize) / binSize) - 1),
		to = binSize * ceiling(max(isize) / binSize),
		by = binSize
	)
	his <- table(cut(isize, breaks))
	
	# Normalize
	his <- his / max(his)
	
	# Trim tail of values < 1% of the max
	rle <- rle(as.logical(his < 0.01))
	if(isTRUE(tail(rle$values, 1))) {
		his <- his[ 1:(length(his) - tail(rle$lengths, 1)) ]
	}
	
	# Bin coordinates
	from <- as.numeric(sub("^\\(([0-9\\.]+),([0-9\\.]+)]$", "\\1", names(his))) + 1L
	to <- as.numeric(sub("^\\(([0-9\\.]+),([0-9\\.]+)]$", "\\2", names(his)))
	
	# Step plot
	x <- as.vector(rbind(from, to))
	y <- rep(as.vector(his), each=2)
}

# YAML graph
lines <- c(
	"id: 'Insert_size'",
	"section_name: 'Insert size'",
	sprintf("description: 'histogram (10-bp bins), computed from the first %g \"proper\" read pairs in the BAM file. In a paired-end sequencing experiment, it corresponds to the estimated size of the (c)DNA fragment whose ends were sequenced, based on the aligned locations of the ends. The tail of bins < 1%% of the most represented size bin is trimmed.'", nReads),
	"plot_type: 'linegraph'",
	"pconfig:",
	"    id: 'Insert_size_linegraph'",
	"    title: 'Insert size'",
	"    xlab: 'Insert size (bp)'",
	"    ylab: 'Frequency (normalized to the most represented size bin)'",
	"    ymin: 0",
	"    ymax: 1",
	"data:",
	sprintf("    '%s': { %s }", sample, paste(sprintf("%i: %g", x, y), collapse=", "))
)
cat(lines, sep="\n")

