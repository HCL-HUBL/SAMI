#!/usr/bin/env Rscript

### Estimates the distribution of insert sizes from the first read pairs in a BAM file (use a **transcriptome-aligned** BAM file for RNA-seq !)
### Prints to stdout the YAML formatted results, for handling by MultiQC (automatically if stored in a file with name "*_mqc.yaml")
### Requires SAMtools to be available from the PATH
### 
### EXAMPLE : ./insertSize.R sampleA ../STAR/sampleA/Aligned.toTranscriptome.out.bam > sampleA_mqc.yaml
### 
### Author : <sylvain.mareschal@lysarc.org>

nReads <- 1e6
nPoints <- 100L



# CLI arguments
args <- commandArgs(TRUE)
if(!length(args) %in% 2:3) stop("USAGE: ./insertSize.R SAMPLE_ID BAM_FILE [ SAMTOOLS ]")
sample <- args[1]
bamFile <- args[2]
samtools <- args[3]

# Check arguments
if(samtools == "" || is.na(samtools)) samtools <- "samtools"
if(!file.exists(bamFile)) stop("BAM_FILE must exist")

# Get ISIZE field for nth first alignments (R2 & proper-paired only)
isize <- system(sprintf("\"%s\" view -f 0x2 -f 0x80 \"%s\" | cut -f9 | head -%i", samtools, bamFile, nReads), intern=TRUE)
isize <- abs(as.integer(isize))

# Moving average <https://stackoverflow.com/a/4862334>
ma <- function(x, n=5) filter(x, rep(1/n, n), sides=2)

# Smoothed histogram
h <- tabulate(isize)
hs <- ma(h, 60)

# Mode
xmod <- which.max(h)
ymod <- max(h)

# 5th percentiles
i <- hs
i[ is.na(i) ] <- 0L
pct <- c(
	head(which(cumsum(i)/sum(i) > 0.05), 1),
	tail(which(rev(cumsum(rev(i))/sum(i) > 0.05)), 1)
)

# Trim lowly represented insert sizes
idx <- 1:length(hs)
rng <- c(
	tail(which(idx < xmod & hs < ymod*0.01), 1),
	head(which(idx > xmod & hs < ymod*0.01)[1], 1)
)

# Decrease graph sampling
x <- round(seq(from=rng[1], to=rng[2], length=nPoints))
y <- signif(1000 * hs[x] / nReads, 3)

# YAML graph
lines <- c(
	"id: 'Insert_size'",
	"section_name: 'Insert size'",
	sprintf("description: 'distribution, smoothed from the first %g \"proper\" read pairs in the BAM file. In a paired-end sequencing experiment, it corresponds to the estimated size of the (c)DNA fragment whose ends were sequenced, based on the aligned locations of the ends. The plot is trimmed to sizes represented by at least 1%% of the frequency of the most represented size.'", nReads),
	"plot_type: 'linegraph'",
	"pconfig:",
	"    id: 'Insert_size_linegraph'",
	"    title: 'Insert size'",
	"    xlab: 'Insert size (bp)'",
	"    ylab: 'Frequency (arbitrary unit)'",
	"data:",
	sprintf("    %s: { %s }", sample, paste(sprintf("%i: %g", x, y), collapse=", "))
)
cat(lines, sep="\n")
