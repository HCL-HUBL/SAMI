#!/usr/bin/env Rscript --vanilla

# Get FASTQ sets from Nextflow
R1 <- strsplit("!{R1.join("|")}", split="|", fixed=TRUE)[[1]]
R2 <- strsplit("!{R2.join("|")}", split="|", fixed=TRUE)[[1]]

# Check type consistency
type <- strsplit("!{type.join("|")}", split="|", fixed=TRUE)[[1]]
type <- unique(type)
if(length(type) != 1L) stop("A sample can't mix single and paired-end FASTQ")
pairedEnd <- type == "paired"

# For each R1/R2 pair
RG <- NULL
for(i in 1:length(R1)) {
	# RG definition for BAM header (current pair)
	x <- sprintf("ID:%s_%i", "!{sample}", i)
	if("!{CN}" != "") x <- c(x, "CN:!{CN}")
	if("!{PL}" != "") x <- c(x, "PL:!{PL}")
	if("!{PM}" != "") x <- c(x, "PM:!{PM}")
	x <- c(x, "SM:!{sample}")

	# Merge with all read pairs
	RG <- c(RG, paste(x, collapse="\t"))
}

# Print final RG to stdout
cat(paste(RG, collapse=" , "))
