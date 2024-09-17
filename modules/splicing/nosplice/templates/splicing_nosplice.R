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
	mtx <- as.matrix(tab[,-(1:2)])
	rownames(mtx) <- paste(sub("^chr", "", tab$"#CHROM"), tab$"POS", sep=":")
	
	return(mtx)
}

# Add nosplice events to the 'events' table
updateEvents <- function(events, splicing, sites) {
	# Intron retention on the left (on known sites and classical splicing only)
	left.nosplice <- with(
		events[ events$class %in% c("annotated", "anchored-left", "plausible") & rownames(events) %in% splicing ,],
		data.frame(
			name = sprintf("%s?:%i-nosplice", left.chrom, left.pos),
			left.chrom = left.chrom,
			left.strand = "?",
			left.pos = left.pos,
			right.chrom = NA,
			right.strand = NA,
			right.pos = NA,
			class = "nosplice-left"
		)
	)
	left.nosplice <- unique(left.nosplice)
	rownames(left.nosplice) <- left.nosplice$name
	left.nosplice$name <- NULL

	# Intron retention on the right (on known sites and classical splicing only)
	right.nosplice <- with(
		events[ events$class %in% c("annotated", "anchored-right", "plausible") & rownames(events) %in% splicing ,],
		data.frame(
			name = sprintf("nosplice-%s?:%i", right.chrom, right.pos),
			left.chrom = NA,
			left.strand = NA,
			left.pos = NA,
			right.chrom = right.chrom,
			right.strand = "?",
			right.pos = right.pos,
			class = "nosplice-right"
		)
	)
	right.nosplice <- unique(right.nosplice)
	rownames(right.nosplice) <- right.nosplice$name
	right.nosplice$name <- NULL

	# Annotated left nosplice
	site <- sprintf("%s:%i", left.nosplice$left.chrom, left.nosplice$left.pos)
	exons <- strsplit(sites[site,"exons.all"], split=",", fixed=TRUE)
	known <- sapply(exons, function(x) { any(grepl("^[0-9]+$", x)) })
	left.nosplice[ known , "class" ] <- "annotated"
	
	# Annotated right nosplice
	site <- sprintf("%s:%i", right.nosplice$right.chrom, right.nosplice$right.pos)
	exons <- strsplit(sites[site,"exons.all"], split=",", fixed=TRUE)
	known <- sapply(exons, function(x) { any(grepl("^[0-9]+$", x)) })
	right.nosplice[ known , "class" ] <- "annotated"
	
	# Add events
	events <- rbind(events, left.nosplice, right.nosplice)
	
	return(events)
}

# Add nosplice events to the 'groups' table
updateGroups <- function(groups, events) {
	# Left splicing sites
	left <- grep("-nosplice$", rownames(events), value=TRUE)
	left.groups <- data.frame(
		site  = sub("?", "", fixed=TRUE, sub("-nosplice$", "", left)),
		event = left,
		side  = "left"
	)
	
	# Right splicing sites
	right <- grep("^nosplice-", rownames(events), value=TRUE)
	right.groups <- data.frame(
		site  = sub("?", "", fixed=TRUE, sub("^nosplice-", "", right)),
		event = right,
		side  = "right"
	)
	
	# Add groups
	groups <- rbind(groups, left.groups, right.groups)
	
	return(groups)
}

# Add sequencing depth in intron as 'I' for 'nosplice' events
updateI <- function(groups, I) {
	
	### LEFT ###
	
	# nosplice events on a left site
	left <- grep("-nosplice$", groups$event)
	
	# Genomic point to use as I
	split <- strsplit(groups[ left , "site" ], split=":")
	chrom <- sapply(split, "[", 1L)
	pos <- as.integer(sapply(split, "[", 2L))
	at <- paste(chrom, pos+2L, sep=":")
	
	# Extra I rows
	left.depth <- depth[ at, colnames(I) ]
	rownames(left.depth) <- groups[ left , "event" ]
	
	
	### RIGHT ###
	
	# nosplice events on a right site
	right <- grep("^nosplice-", groups$event)
	
	# Genomic point to use as I
	split <- strsplit(groups[ right , "site" ], split=":")
	chrom <- sapply(split, "[", 1L)
	pos <- as.integer(sapply(split, "[", 2L))
	at <- paste(chrom, pos-2L, sep=":")
	
	# Extra I rows
	right.depth <- depth[ at, colnames(I) ]
	rownames(right.depth) <- groups[ right , "event" ]
	
	
	### MERGE ###
	
	# Add I
	I <- rbind(I, left.depth, right.depth)
	
	# Check consistency
	if(!identical(rownames(I), groups$event)) stop("Inconsistency")
	
	return(I)
}

# Recompute 'S' from scratch, including 'nosplice' counts
updateS <- function(groups, I, S) {
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
	sums <- sums[ paste(groups$site, groups$side) ,]
	rownames(sums) <- groups$event
	
	# Check
	if(!identical(dimnames(sums), dimnames(I))) stop("Inconsistency")
	
	# Compute 'S' from 'I' and totals
	S <- sums - I
	
	return(S)
}



timedMessage("Parsing RDS files...")

I <- readRDS("I.rds")
S <- readRDS("S.rds")
groups <- readRDS("groups.rds")
sites <- readRDS("sites.rds")
events <- readRDS("events.rds")
splicing <- readRDS("splicing.rds")

timedMessage("Parsing depth files...")

depth <- parseDepth("depth.bed")

timedMessage("Adding events...")

events <- updateEvents(events, splicing, sites)

timedMessage("Adding groups...")

groups <- updateGroups(groups, events)

timedMessage("Adding I...")

I <- updateI(groups, I)

timedMessage("Recomputing S...")

S <- updateS(groups, I, S)

timedMessage("Resorting...")

o <- order(groups$site, groups$event, groups$side)
groups <- groups[ o ,]
rownames(groups) <- NULL
I <- I[ o ,]
S <- S[ o, ]

timedMessage("Exporting...")

dir.create("out")
saveRDS(I, file="out/I.rds")
saveRDS(S, file="out/S.rds")
saveRDS(groups, file="out/groups.rds")
saveRDS(sites, file="out/sites.rds")
saveRDS(events, file="out/events.rds")
saveRDS(depth, file="out/depth.rds")

timedMessage("done")
