#!/usr/bin/env Rscript --vanilla

# Collect Nextflow arguments
ncores <- as.integer("!{task.cpus}")
geneFile <- "!{genes}"     
exonFile <- "!{exons}"
intronFile <- "!{introns}"
chromosomes <- strsplit("!{chromosomes}", split=",")[[1]]
min.reads.unknown <- as.integer("!{min_reads_unknown}")
transcriptFile <- "transcripts.tsv"
stranded <- "!{stranded}"



# Print a log message with date and time
timedMessage <- function(...) {
	message(Sys.time(), " : ", ...)
}

# Collect STAR junctions from individual RDT files into a single filtered matrix : mtx[ event , sample ] = read.count
collectJunctions <- function(chromosomes, stranded) {
	# Storage
	lst <- list()

	samples <- NULL
	# Regular junctions
	files <- dir("gapFiles", pattern="\\.DNA\\.MD\\.sort\\.bam\\.tsv$", full.names=TRUE)
	for(file in files) {
		# Parse file
		message("- ", file)
		tab <- read.table(
			file, sep="\t", quote=NULL, comment.char="",
			col.names = c("chrom", "start", "end", "bases", "orientation", "reads"),
			colClasses = c("character", "integer", "integer", "character", "character", "integer"),
		)
		
		# Strand
		tab$strand <- "?"
		
		# Trust most common splicing sites
		tab[ tab$bases == "CT-AC" , "strand" ] <- "-"
		tab[ tab$bases == "GT-AG" , "strand" ] <- "+"
		
		# But trust more correct pair orientation, if available
		if(stranded == "R1") {
			tab[ tab$orientation == "F1R2" , "strand" ] <- "+"
			tab[ tab$orientation == "F2R1" , "strand" ] <- "-"
		} else if(stranded == "R2") {
			tab[ tab$orientation == "F1R2" , "strand" ] <- "-"
			tab[ tab$orientation == "F2R1" , "strand" ] <- "+"
		}
		
		# Reshape and filter chromosome
		tab$chrom <- factor(sub("^chr", "", tab$chrom), levels=chromosomes)
		tab <- tab[ !is.na(tab$chrom) ,]
		
		# Genomically ordered ID (SJ respects that, chimeric doesn't)
		tab$left <- pmin(tab$start, tab$end)
		tab$right <- pmax(tab$start, tab$end)
		tab$ID <- with(tab, sprintf("%s%s:%i-%s%s:%i", chrom, strand, left, chrom, strand, right))
		
		# Gather
		sample <- sub("\\.DNA\\.MD\\.sort\\.bam\\.tsv$", "", basename(file))
		samples <- c(samples, sample)
		tab$sample <- sample
		lst[[ file ]] <- tab[, c("sample", "ID", "reads") ]
	}
	
	# Chimeric junctions
	files <- dir("chimericFiles", pattern="_Chimeric\\.out\\.junction$", full.names=TRUE)
	for(file in files) {
		# Parse file
		message("- ", file)
		chi <- read.table(
			file, sep="\t", quote=NULL, comment.char="", skip=1,
			col.names  = c("A.chrom",   "A.break", "A.strand",  "B.chrom",   "B.break", "B.strand",  "type",    "A.rep",   "B.rep",   "read",      "A.start", "A.CIGAR",   "B.start", "B.CIGAR",   "num_chim_aln", "max_poss_aln_score", "non_chim_aln_score", "this_chim_aln_score", "bestall_chim_aln_score", "PEmerged_bool", "RG"),
			colClasses = c("character", "integer", "character", "character", "integer", "character", "integer", "integer", "integer", "character", "integer", "character", "integer", "character", "integer",      "integer",            "integer",            "integer",             "integer",                "integer",       "character"),
			fill=TRUE
		)
		
		# Recurrence
		chi <- chi[,1:6]
		id <- apply(chi, 1, paste, collapse="|")
		i <- !duplicated(id)
		chi <- chi[i,]
		chi$reads <- as.integer(table(id)[ id[i] ])
		
		# Reshape and filter chromosome
		chi$A.chrom <- factor(sub("^chr", "", chi$A.chrom), levels=chromosomes)
		chi$B.chrom <- factor(sub("^chr", "", chi$B.chrom), levels=chromosomes)
		chi <- chi[ !is.na(chi$A.chrom) & !is.na(chi$B.chrom) ,]
		
		# Genomically ordered ID
		A.left <- rep(NA, nrow(chi))
		for(i in 1:nrow(chi)) A.left[i] <- as.integer(chi$A.chrom[i]) <= as.integer(chi$B.chrom[i]) || chi$A.break[i] <= chi$B.break[i]
		chi$ID <- ""
		chi$ID[  A.left ] <- with(chi[  A.left ,], sprintf("%s%s:%i-%s%s:%i", A.chrom, A.strand, A.break, B.chrom, B.strand, B.break))
		chi$ID[ !A.left ] <- with(chi[ !A.left ,], sprintf("%s%s:%i-%s%s:%i", A.chrom, A.strand, A.break, B.chrom, B.strand, B.break))
		
		# Gather
		sample <- sub("_Chimeric\\.out\\.junction$", "", basename(file))
		samples <- c(sample, samples)
		chi$sample <- sample
		lst[[ file ]] <- chi[, c("sample", "ID", "reads") ]
	}

	message("- Reshaping")

	# Prepare a storage matrix
	IDs <- sort(unique(unlist(lapply(lst, "[[", "ID"))))
	samples <- sort(unique(samples))
	mtx <- matrix(
		data = 0L,
		nrow = length(IDs),
		ncol = length(samples),
		dimnames = list(IDs, samples)
	)
	
	# Transfer from the list to the matrix (additive)
	for(chunk in lst) {
		mtx[ as.matrix(chunk[,c("ID", "sample")]) ] <- mtx[ as.matrix(chunk[c("ID","sample")]) ] + chunk$reads
	}
	
	return(mtx)
}

# Group splicing event IDs (chrX:NNN-NNN) which share a common splicing site
sharedSites <- function(events) {
	# Splicing sites (either genomic left or genomic right)
	left <- paste(events$left.chrom, events$left.pos, sep=":")
	right <- paste(events$right.chrom, events$right.pos, sep=":")
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
annotateSingleSite <- function(site, events.indexes, mtx, events, genes, exons, preferred) {
	
	# Determine whether the considered genomic position corresponds to an exon boundary or not
	annotateExon <- function(chrom, position, genes, exons, preferred) {
		anno <- NULL
		
		# Exons overlapping the position of interest
		overlap <- exons$slice(chrom, position - 1L, position + 1L)
		
		# Position of interest is exon start
		ovl <- overlap[ overlap$start - 1L == position ,]
		if(nrow(ovl) > 0L) {
			ovl <- ovl[, c("transcript", "symbol", "strand", "groupPosition", "groupPosition") ]
			colnames(ovl) <- c("transcript", "gene", "strand", "exon", "anno")
			ovl$transcript <- sub("\\..+$", "", ovl$transcript)
			ovl$anno <- sprintf("[%i", ovl$anno)
			anno <- rbind(anno, ovl)
		}
		
		# Position of interest is exon end
		ovl <- overlap[ overlap$end + 1L == position ,]
		if(nrow(ovl) > 0L) {
			ovl <- ovl[, c("transcript", "symbol", "strand", "groupPosition", "groupPosition") ]
			colnames(ovl) <- c("transcript", "gene", "strand", "exon", "anno")
			ovl$transcript <- sub("\\..+$", "", ovl$transcript)
			ovl$anno <- sprintf("%i]", ovl$anno)
			anno <- rbind(anno, ovl)
		}
		
		# Position of interest is inside exon
		ovl <- overlap[ overlap$start - 1L != position & overlap$end + 1L != position ,]
		if(nrow(ovl) > 0L) {
			ovl <- ovl[, c("transcript", "symbol", "strand", "groupPosition", "groupPosition") ]
			colnames(ovl) <- c("transcript", "gene", "strand", "exon", "anno")
			ovl$transcript <- sub("\\..+$", "", ovl$transcript)
			anno <- rbind(anno, ovl)
		}
		
		# Merge
		if(is.null(anno)) {
			exon.all <- NA
		} else {
			exon.all <- paste(unique(anno[ order(anno$exon) , "anno" ]), collapse=",")
		}
		
		# Preferred transcript
		trans.prf <- intersect(preferred, anno$transcript)
		if(length(trans.prf) > 0L) {
			if(is.null(anno)) {
				exon.prf <- NA
			} else {
				exon.prf <- paste(anno[ match(trans.prf, anno$transcript) , "anno" ], collapse=",")
			}
		} else {
			exon.prf <- NA
		}
		
		# Transcript list
		transcripts <- sort(unique(anno$transcript))
		
		# No exon : looking for an intron or UTR
		if(is.null(anno)) {
			anno <- genes$slice(chrom, position - 1L, position + 1L)
			column <- "name"
		} else {
			column <- "gene"
		}
		
		# Gene list (with strand)
		site.genes <- sort(unique(sprintf("%s (%s)", anno[[column]], anno$strand)))
		
		return(
			list(
				exon.all = exon.all,
				exon.preferred = exon.prf,
				transcripts = paste(transcripts, collapse=","),
				transcripts.preferred = paste(trans.prf, collapse=","),
				genes = paste(site.genes, collapse=",")
			)
		)
	}
	
	# Reads supporting inclusion
	I <- mtx[ events.indexes , , drop=FALSE ]
	
	# Reads supporting exclusion
	S <- t(apply(I, 2, sum) - t(I))
	
	# Events involving the site of interest
	groups <- data.frame(
		site = site,
		event = rownames(events)[ events.indexes ],
		stringsAsFactors = FALSE
	)
	
	# Identify left and rigth sites of events
	groups[ groups$site == sub("^([^:]+)[+?-]:([0-9]+)-([^:]+)[+?-]:([0-9]+)$", "\\1:\\2", groups$event) , "side" ] <- "left"
	groups[ groups$site == sub("^([^:]+)[+?-]:([0-9]+)-([^:]+)[+?-]:([0-9]+)$", "\\3:\\4", groups$event) , "side" ] <- "right"

	# Site of interest
	chrom <- sub("^([^:]+):([0-9]+)$", "\\1", site)
	position <- as.integer(sub("^([^:]+):([0-9]+)$", "\\2", site))
	site.exons <- annotateExon(chrom, position, genes, exons, preferred)
	sites <- data.frame(
		row.names = site,
		chrom = chrom,
		position = position,
		reads = sum(I),
		genes = site.exons$genes,
		transcripts = site.exons$transcripts,
		exons.all = site.exons$exon.all,
		transcripts.preferred = site.exons$transcripts.preferred,
		exons.preferred = site.exons$exon.preferred,
		stringsAsFactors = FALSE
	)
	
	return(
		list(
			I = I,
			S = S,
			sites = sites,
			groups = groups
		)
	)
}

# Process a subset of the whole sameSite list, and return data.frame chunks as a list (for efficient parallelization)
annotateSiteChunk <- function(sites.events, mtx, events, genes, exons, annotateSingleSite, preferred) {
	# Dependencies
	library(Rgb)
	
	# Process sites sequentially
	out <- vector(mode="list", length=length(sites.events))
	for(i in 1:length(sites.events)) {
		out[[i]] <- annotateSingleSite(names(sites.events)[i], sites.events[[i]], mtx, events, genes, exons, preferred)
	}
	
	# rbind elements separately
	mrg <- list()
	for(element in names(out[[1]])) {
		mrg[[element]] <- do.call(rbind, lapply(out, "[[", element))
	}
	
	return(mrg)
}

# Process all of the sameSite list, using parallelization
annotateAllSites <- function(sites.events, ncores, mtx, events, genes, exons, preferred) {
	# Create a cluster for parallelization
	cluster <- makeCluster(spec=ncores)
	
	message("- Submitting ", length(sites.events), " sites on ", ncores, " CPUs...")
	out <- clusterApply(
		cl = cluster,
		fun = annotateSiteChunk,
		x = split(sites.events, 1:ncores),
		mtx = mtx,
		events = events,
		genes = genes,
		exons = exons,
		preferred = preferred,
		annotateSingleSite = annotateSingleSite
	)
	
	# Close the cluster
	stopCluster(cluster)
	
	# rbind elements separately
	mrg <- list()
	for(element in names(out[[1]])) {
		mrg[[element]] <- do.call(rbind, lapply(out, "[[", element))
	}
	
	return(mrg)
}

# Classify genomic junctions according to known transcripts
classifyJunctions <- function(ID, introns, exons, chromosomes) {
	# Coordinates of events
	regex <- "^([^:]+)([+?-]):([0-9]+)-([^:]+)([+?-]):([0-9]+)$"
	events <- data.frame(row.names=ID, stringsAsFactors=FALSE)
	events$left.chrom   <- factor(sub(regex, "\\1", ID), levels=chromosomes)
	events$left.strand  <- sub(regex, "\\2", ID)
	events$left.pos     <- as.integer(sub(regex, "\\3", ID))
	events$right.chrom  <- factor(sub(regex, "\\4", ID), levels=chromosomes)
	events$right.strand <- sub(regex, "\\5", ID)
	events$right.pos    <- as.integer(sub(regex, "\\6", ID))
	
	# Start considering all junctions as random
	events$class <- "unknown"

	# Check if junction binds two exons
	exonStarts  <- sprintf("%s:%i", exons$extract(,"chrom"), exons$extract(,"start") - 1L)
	exonEnds    <- sprintf("%s:%i", exons$extract(,"chrom"), exons$extract(,"end") + 1L)
	eventLefts  <- sprintf("%s:%i", events$left.chrom, events$left.pos)
	eventRights <- sprintf("%s:%i", events$right.chrom, events$right.pos)
	leftMatch   <- eventLefts %in% c(exonStarts, exonEnds)
	rightMatch  <- eventRights %in% c(exonStarts, exonEnds)
	events[ leftMatch | rightMatch , "class" ] <- "anchored"
	events[ leftMatch & rightMatch , "class" ] <- "plausible"

	# Check if junction exists in annotation
	events[ events$left.chrom == events$right.chrom & sprintf("%s:%i-%i", events$left.chrom, events$left.pos, events$right.pos) %in% introns, "class" ] <- "annotated"
	
	return(events)
}

# Remove poor junctions
filterJunctions <- function(events, mtx, min.reads.unknown=10L) {
	# Apply a minimum read count filter to "unknown" junctions only
	filter.1 <- events$class != "unknown" | apply(mtx, 1, max, na.rm=TRUE) >= min.reads.unknown
	message("- ", sum(!filter.1), " / ", nrow(mtx), " events filtered out (unknown with less than ", min.reads.unknown, " reads in each sample)")
	
	# Discard left = right events
	filter.2 <- events$left.chrom != events$right.chrom | events$left.pos != events$right.pos
	message("- ", sum(!filter.2), " / ", nrow(mtx), " events filtered out (left = right)")
	
	return(filter.1 & filter.2)
}

# Parse a table of preferred transcripts
parsePreferred <- function(file, exons) {
	# Parse
	tmp <- scan(file, what=list("", ""), sep="\t", quiet=TRUE)
	x <- sub("\\.[0-9]+$", "", tmp[[2]])
	names(x) <- tmp[[1]]
	
	# Check existence in annotation
	known <- x %in% sub("\\.[0-9]+$", "", exons$extract(,"transcript"))
	if(any(!known)) stop("Preferred transcripts not found in annotation : ", paste(x[!known], collapse=", "))
	
	return(x)
}



timedMessage("Loading dependencies...")

library(Rgb)
library(parallel)

timedMessage("Parsing annotation...")

genes <- readRDT(geneFile)
exons <- readRDT(exonFile)
introns <- readRDS(intronFile)

timedMessage("Parsing junction files...")

mtx <- collectJunctions(chromosomes, stranded)

timedMessage("Classifying...")

events <- classifyJunctions(rownames(mtx), introns, exons, chromosomes)

timedMessage("Filtering...")

filter <- filterJunctions(events, mtx, min.reads.unknown)
if(!any(filter)) stop("No splicing event to work with")
mtx <- mtx[ filter ,, drop=FALSE ]
events <- events[ filter ,, drop=FALSE ]

timedMessage("Grouping per site...")

sites.events <- sharedSites(events)

timedMessage("Parsing preferred transcript file...")

preferred <- parsePreferred(transcriptFile, exons)

timedMessage("Annotating...")

out <- annotateAllSites(sites.events, ncores, mtx, events, genes, exons, preferred)

timedMessage("Exporting...")

saveRDS(out$I, file="I.rds")
saveRDS(out$S, file="S.rds")
saveRDS(out$groups, file="groups.rds")
saveRDS(out$sites, file="sites.rds")
saveRDS(events, file="events.rds")

timedMessage("done")
