#!/usr/bin/env Rscript --vanilla

# Print a log message with date and time
timedMessage <- function(...) {
	message(Sys.time(), " : ", ...)
}

# Collect depth values as a matrix
parseDepth <- function(bedFile) {
	# Parse depth files
	tab <- read.table(bedFile, sep="\t", header=TRUE, check.names=FALSE, comment.char="")
	
	# Fix sample names
	colnames(tab) <- sub("\\.DNA\\.MD\\.sort\\.bam$", "", basename(colnames(tab)))
	
	# Reformat as a matrix
	mtx <- as.matrix(tab[ , -(1:2) , drop=FALSE ])
	rownames(mtx) <- paste(sub("^chr", "", tab$"#CHROM"), tab$"POS", sep=":")
	
	return(mtx)
}

# Add sequencing depth in intron as 'I' for 'nosplice' events
updateI <- function(groups, I, depth, shift=3L) {
	
	### LEFT ###
	
	# nosplice events on a left site
	left <- grep("-nosplice$", groups$event)
	
	# Genomic point to use as I
	split <- strsplit(groups[ left , "site" ], split=":")
	chrom <- sapply(split, "[", 1L)
	pos <- as.integer(sapply(split, "[", 2L))
	at <- paste(chrom, pos + shift - 1L, sep=":")
	
	# Update I
	I[ groups[ left , "event" ] ,] <- depth[ at, colnames(I) ]
	
	
	### RIGHT ###
	
	# nosplice events on a right site
	right <- grep("^nosplice-", groups$event)
	
	# Genomic point to use as I
	split <- strsplit(groups[ right , "site" ], split=":")
	chrom <- sapply(split, "[", 1L)
	pos <- as.integer(sapply(split, "[", 2L))
	at <- paste(chrom, pos - shift + 1L, sep=":")
	
	# Update I
	I[ groups[ right , "event" ] ,] <- depth[ at, colnames(I) ]
	
	
	return(I)
}

# Compute 'S' from scratch, including 'nosplice' counts
computeS <- function(groups, I) {
	# Cluster I values by splicing site
	clusters <- tapply(INDEX=paste(groups$site, groups$side), X=1:nrow(groups), FUN=c)
	
	# Storage for sums
	sums <- matrix(
		as.integer(NA),
		nrow = length(clusters),
		ncol = ncol(I),
		dimnames = list(
			names(clusters),
			colnames(I)
		)
	)
	
	# Compute sums
	for(i in 1:nrow(sums)) sums[i,] <- apply(I[ clusters[[i]] ,, drop=FALSE ], 2, sum)
	
	# Reshape as 'S'
	sums <- sums[ paste(groups$site, groups$side) ,, drop=FALSE ]
	rownames(sums) <- groups$event
	
	# Check
	if(!identical(dimnames(sums), dimnames(I))) stop("Inconsistency")
	
	# Compute 'S' from 'I' and totals
	S <- sums - I
	
	return(S)
}



timedMessage("Parsing RDS files...")

I <- readRDS("I.rds")
groups <- readRDS("groups.rds")
sites <- readRDS("sites.rds")
events <- readRDS("events.rds")

timedMessage("Parsing depth files...")

depth <- parseDepth("depth.bed")

timedMessage("Updating I...")

I <- updateI(groups, I, depth)

timedMessage("Computing S...")

S <- computeS(groups, I)

timedMessage("Exporting...")

dir.create("out")
saveRDS(I, file="out/I.rds")
saveRDS(S, file="out/S.rds")
saveRDS(groups, file="out/groups.rds")
saveRDS(sites, file="out/sites.rds")
saveRDS(events, file="out/events.rds")

timedMessage("done")
