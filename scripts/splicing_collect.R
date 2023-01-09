#!/usr/bin/env Rscript

# Collect CLI arguments
args <- commandArgs(TRUE)
if(length(args) < 7L) stop("USAGE : ./splicing_collect.R NCORES exons.rdt introns.rds output.rds CHROMOSOMES junctions_1.rdt [ junctions_2.rdt [ ... ] ]")
ncores <- as.integer(args[1])
exonFile <- args[2]
intronFile <- args[3]
outputFile <- args[4]
chromosomes <- strsplit(args[5], split=",")[[1]]
junctionFiles <- args[ 6:length(args) ]



# Collect STAR junctions from individual RDT files into a single filtered matrix : mtx[ event , sample ] = read.count
collectJunctions <- function(files) {
	# Collect junctions as a list
	lst <- list()
	for(file in files) {
		# Extract junctions
		tmp <- readRDT(file)
		tab <- tmp$extract(, c("chrom", "start", "end", "reads"))
		
		# Genomically ordered ID (SJ respects that, chimeric doesn't)
		tab$left <- pmin(tab$start, tab$end)
		tab$right <- pmax(tab$start, tab$end)
		tab$ID <- with(tab, sprintf("chr%s:%i-%i", chrom, left, right))
		
		# Gather
		lst[[ tmp$name ]] <- tab[, c("ID", "reads") ]
	}
	
	# Prepare a storage matrix
	IDs <- unique(unlist(lapply(lst, "[[", "ID")))
	mtx <- matrix(
		data = 0L,
		nrow = length(IDs),
		ncol = length(lst),
		dimnames = list(IDs, names(lst))
	)
	
	# Transfer from the list to the matrix
	for(sample in names(lst)) {
		mtx[ lst[[sample]]$ID , sample ] <- lst[[sample]]$reads
	}
	
	return(mtx)
}

# Group splicing event IDs (chrX:NNN-NNN) which share a common splicing site
sharedSites <- function(IDs) {
	# Splicing sites (either genomic left or genomic right)
	left <- sub("^chr([0-9XY]+):([0-9]+)-([0-9]+)$", "\\1:\\2", IDs)
	right <- sub("^chr([0-9XY]+):([0-9]+)-([0-9]+)$", "\\1:\\3", IDs)
	sites <- unique(c(left, right))
	
	# For each junction, get the indexes of left and right sites in the site dictionnary
	left.site <- match(left, sites)
	right.site <- match(right, sites)
	
	# For each site in the site dictionnary, get the indexes of junctions whose left or right boundary matches
	site.lefts  <- tapply(X=1:length(left),  INDEX=left.site,  FUN=c)
	site.rights <- tapply(X=1:length(right), INDEX=right.site, FUN=c)
	
	# Concatenate left and right matches into a single list
	site.lefts  <- site.lefts [ as.character(1:length(sites)) ]
	site.rights <- site.rights[ as.character(1:length(sites)) ]
	sameSite <- mapply(FUN=c, site.lefts, site.rights)
	names(sameSite) <- sites
	
	return(sameSite)
}

# Aggregate information as a data.frame about multiple events involved in a single site
annotateSingleSite <- function(site, events.indexes, mtx, exons, readThrough=FALSE, chromosomes=c(1:22, "X", "Y")) {
	# Reads supporting inclusion
	I <- mtx[ events.indexes , , drop=FALSE ]
	
	# Position of the splicing event
	chrom <- sub("^([0-9XY]+):([0-9]+)$", "\\1", site)
	
	# Reads supporting exclusion
	S <- t(apply(I, 2, sum) - t(I))
	PSI <- I/(I+S)
	
	# Annotation of events involving the site
	IDs <- rownames(mtx)[ events.indexes ]
	chrom <- sub("^(chr)?([^:]+):([0-9]+)-([0-9]+)$", "\\2", IDs)
	left <- as.integer(sub("^(chr)?([^:]+):([0-9]+)-([0-9]+)$", "\\3", IDs))
	right <- as.integer(sub("^(chr)?([^:]+):([0-9]+)-([0-9]+)$", "\\4", IDs))
	if(any(left > right)) stop("left / right inconsistency")
	
	# Liste all genes involved
	symbol <- list()
	transcript <- list()
	isReadThrough <- NULL
	for(j in 1:length(events.indexes)) {
		# Genes at left position
		g <- exons$slice(chrom=chrom[j], start=left[j]-10L, end=left[j]+10L)
		symbol.left <- unique(sprintf("%s (%s)", g$symbol, g$strand))
		transcript.left <- unique(g$transcript)
		
		# Genes at right position
		g <- exons$slice(chrom=chrom[j], start=right[j]-10L, end=right[j]+10L)
		symbol.right <- unique(sprintf("%s (%s)", g$symbol, g$strand))
		transcript.right <- unique(g$transcript)
		
		# Read-through (2 distinct genes at splicing start and end)
		isReadThrough[j] <- length(symbol.left) > 0L & length(symbol.right) > 0L & length(intersect(symbol.left, symbol.right)) == 0L
		
		# Pooled genes
		symbol[[j]] <- union(symbol.left, symbol.right)
		transcript[[j]] <- union(transcript.left, transcript.right)
	}
	
	# Output chunk (one row per event)
	tab <- data.frame(
		site = site,
		site.reads = sum(I),
		ID = sub("chr", "", rownames(I)),
		ID.symbol = sapply(symbol, paste, collapse=", "),
		ID.transcript = sapply(transcript, paste, collapse=", "),
		readThrough = isReadThrough,
		I = I,
		S = S,
		PSI = PSI,
		depth = I+S,
		max.PSI = apply(PSI, 1, max, na.rm=TRUE),
		chrom = factor(chrom, levels=chromosomes),
		left = left,
		right = right,
		stringsAsFactors = FALSE,
		row.names = NULL,
		check.names = FALSE
	)
	
	# Discard read-through alternatives
	if(!isTRUE(readThrough)) tab <- tab[ !tab$readThrough ,]
	
	# Discard left = right events
	tab <- tab[ tab$left != tab$right ,]
	
	# Additional data
	tab$site.is.ID <- ifelse(tab$left == as.integer(sub("^([0-9XY]+):([0-9]+)$", "\\2", tab$site)), "left", "right")
	
	return(tab)
}

# Process a subset of the whole sameSite list, and return data.frame chunks as a list (for efficient parallelization)
annotateSiteChunk <- function(sameSite, mtx, exons, annotateSingleSite, ...) {
	# Dependencies
	library(Rgb)
	
	# Process sites sequentially
	out <- vector(mode="list", length=length(sameSite))
	for(i in 1:length(sameSite)) {
		out[[i]] <- annotateSingleSite(names(sameSite)[i], sameSite[[i]], mtx, exons, ...)
	}
	
	# Merge rows
	tab <- do.call(rbind, out)
	
	return(tab)
}

# Process all of the sameSite list, using parallelization
annotateAllSites <- function(sameSite, ncores, mtx, exons, ...) {
	# Create a cluster for parallelization
	cluster <- makeCluster(spec=ncores)
	
	message("- Submitting ", length(sameSite), " sites on ", ncores, " CPUs...")
	out <- clusterApply(
		cl = cluster,
		fun = annotateSiteChunk,
		x = split(sameSite, 1:ncores),
		mtx = mtx,
		exons = exons,
		annotateSingleSite = annotateSingleSite,
		...
	)
	
	# Close the cluster
	stopCluster(cluster)

	message("- Merging...")
	tab <- do.call(rbind, out)
	
	return(tab)
}

# Classify genomic junctions according to known transcripts
classifyJunctions <- function(events, introns, exons) {
	# Start considering all junctions as random
	events$class <- "unknown"

	# Check if junction binds two exons
	exonStarts  <- sprintf("%s:%i", exons$extract(,"chrom"), exons$extract(,"start") - 1L)
	exonEnds    <- sprintf("%s:%i", exons$extract(,"chrom"), exons$extract(,"end") + 1L)
	eventLefts  <- sprintf("%s:%i", events$chrom, events$left)
	eventRights <- sprintf("%s:%i", events$chrom, events$right)
	leftMatch   <- eventLefts %in% c(exonStarts, exonEnds)
	rightMatch  <- eventRights %in% c(exonStarts, exonEnds)
	events[ leftMatch | rightMatch , "class" ] <- "anchored"
	events[ leftMatch & rightMatch , "class" ] <- "plausible"

	# Check if junction exists in annotation
	events[ events$ID %in% introns, "class" ] <- "annotated"
	
	return(events$class)
}



message("Loading dependencies...")

library(Rgb)
library(parallel)

message("Parsing annotation...")

exons <- readRDT(exonFile)
introns <- readRDS(intronFile)

message("Parsing junction files...")

mtx <- collectJunctions(junctionFiles)

message("Grouping per site...")

sameSite <- sharedSites(rownames(mtx))

message("Annotating...")

events <- annotateAllSites(sameSite, ncores, mtx, exons, chromosomes=chromosomes)

message("Classifying...")

events$class <- classifyJunctions(events, introns, exons)

message("Exporting...")

saveRDS(events, file=outputFile)

message("done")
