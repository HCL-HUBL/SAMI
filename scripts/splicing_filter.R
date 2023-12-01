#!/usr/bin/env Rscript --vanilla

# Collect CLI arguments
args <- commandArgs(TRUE)
if(length(args) != 10L) stop("USAGE : ./splicing_filter.R NCORES target.gtf exons.rdt PLOT FUSIONS MIN_I MIN_PSI SYMBOLS|all CLASSES FOCUS")
ncores <- as.integer(args[1])
targetFile <- args[2]
exonFile <- args[3]
plot <- as.logical(args[4])
fusions <- as.logical(args[5])
min.I <- as.integer(args[6])
min.PSI <- as.double(args[7])
symbols <- args[8]
classes <- args[9]
classList <- strsplit(classes, split=",")[[1]]
focus <- args[10]
if(identical(focus, "none")) { focusList <- NULL
} else                       { focusList <- strsplit(focus, split=",")[[1]]
}



# Print a log message with date and time
timedMessage <- function(...) {
	message(Sys.time(), " : ", ...)
}

# Filter events of classes of interest
filterClass <- function(events, classes=classList) {
	# Already computed, just to filter
	events$filter.class <- events$class %in% classes
	
	return(events)
}

# Filter events with any site overlapping a gene of interest
filterSymbol <- function(events, groups, sites, symbols=NULL) {
	# Site level
	if(length(symbols) > 0L) {
		sites.genes <- strsplit(gsub(" \\([+-]\\)", "", sites$genes), split=",")
		sites$filter.symbol <- sapply(sites.genes, function(x) { any(x %in% symbols) })
	} else {
		sites$filter.symbol <- TRUE
	}
	
	# Group level
	mrg <- merge(x=groups, y=sites[,"filter.symbol",drop=FALSE], by.x="site", by.y="row.names", all=TRUE, sort=FALSE)
	if(!identical(mrg[,1:3], groups[,1:3])) stop("Merging error")
	groups <- mrg
	
	# Event level (left OR right)
	events.filter.symbol <- tapply(X=groups$filter.symbol, INDEX=groups$event, FUN=any)
	events[ , "filter.symbol" ] <- as.logical(events.filter.symbol[ rownames(events) ])

	return(list(events=events, groups=groups, sites=sites))
}

# Filter events without enough supporting reads
filterI <- function(events, I, min.I) {
	# Deduplicate I at 'event' x 'sample' level
	events.I <- I[ rownames(events) ,]

	# Events with at least min.I supporting reads for each sample
	events.filter.I <- events.I >= min.I
	
	return(events.filter.I)
}

# Filter events representing a minor alternative at any splicing sites
filterPSI <- function(events, groups, I, S, min.PSI) {
	# PSI at 'groups' x 'sample' level
	PSI <- I / (I + S)

	# PSI filtering at 'groups' x 'sample' level
	groups.filter.PSI <- !is.na(PSI) & PSI >= min.PSI

	# PSI filtering at 'event' x 'sample' level (left AND right)
	left  <- groups.filter.PSI[ groups$side == "left" ,, drop=FALSE ][ rownames(events) ,, drop=FALSE ]
	right <- groups.filter.PSI[ groups$side == "right" ,, drop=FALSE ][ rownames(events) ,, drop=FALSE ]
	events.filter.PSI <- left & right
	
	return(events.filter.PSI)
}

# Pre-filter data to minimize transfers during parallelization
preparePlots <- function(candidates, events, groups, sites, events.filter.all) {
	# All genes involved for each candidate
	left.genes <- strsplit(candidates$left.genes, ",")
	right.genes <- strsplit(candidates$right.genes, ",")
	genes <- mapply(left.genes, right.genes, FUN=function(x, y) unique(c(x, y)))

	# All plots to produce
	n <- sapply(genes, length)
	toPlot <- data.frame(
		symbol = sub(" \\([+-]\\)$", "", unlist(genes)),
		sample = rep(candidates$sample, n)
	)
	
	# Reduce redundancy
	toPlot <- unique(toPlot)
	toPlot <- as.list(toPlot)
	
	message("- ", length(toPlot$symbol), " plots across ", length(unique(toPlot$symbol)), " genes and ", length(unique(toPlot$sample)), " samples")
	
	# Pre-filter exons
	toPlot$exons <- list()
	for(symbol in unique(toPlot$symbol)) {
		gene <- exons$extract(exons$extract(,"symbol") == symbol)
		for(i in which(toPlot$symbol == symbol)) toPlot$exons[[i]] <- gene
	}
	
	# Pre-filter events based on symbol
	toPlot$events <- list()
	for(symbol in unique(toPlot$symbol)) {
		# Events macthing symbol on either site
		match.symbol <- sapply(
			strsplit(gsub(" \\([+-]\\)", "", sites$genes), split=", "),
			`%in%`, x=symbol
		)
		evt <- events[ unique(groups[ groups$site %in% rownames(sites)[ match.symbol ] , "event" ]) ,]
		
		# Genes (symbol only)
		evt$left.genes <- gsub(" \\([+-]\\)", "", sites[ sprintf("%s:%i", evt$left.chrom, evt$left.pos) , "genes" ])
		evt$right.genes <- gsub(" \\([+-]\\)", "", sites[ sprintf("%s:%i", evt$right.chrom, evt$right.pos) , "genes" ])
		
		# Storage
		for(i in which(toPlot$symbol == symbol)) toPlot$events[[i]] <- evt
	}
	
	# Pre-filter events based on support in sample
	for(i in 1:length(toPlot$events)) {
		# Mark alternatives passing all filters for the considered sample
		eventNames <- rownames(toPlot$events[[i]])
		sample <- toPlot$sample[i]
		toPlot$events[[i]]$filter <- events.filter.all[ eventNames , sample ]
		
		# Supporting reads
		toPlot$events[[i]]$reads <- I[ eventNames , sample ]
		
		# Restrict to events supported in the considered sample
		toPlot$events[[i]] <- toPlot$events[[i]][ toPlot$events[[i]]$reads > 0L ,]
		
		# Get candidates' IDs
		toPlot$events[[i]] <- merge(x=toPlot$events[[i]], y=candidates[ candidates$sample == sample , c("junction", "ID") ], by.x="row.names", by.y="junction", all.x=TRUE, all.y=FALSE)
	}
	
	return(toPlot)
}

# Plots junctions over a simplified representation of the transcript (normalized exon and intron sizes)
plot.normalized <- function(evt, sample, symbol, exons, outDir="out", bamDir="out/BAM", trackDir="out/depth", shape=1) {
	
	library(Rgb)
	
	# Produce (or retrieve existing) track.table objects with sequencing depth in a specific locus
	depth <- function(sample, bamFile, chrom, start, end, trackDir="out/depth", qBase=30, qMap=30, chromosomes=c(1:22, "X", "Y"), assembly="GRCh38") {
		# Depth binary
		trackFile <- sprintf("%s/%s_%s-%i-%i.rdt", trackDir, sample, chrom, start, end)
		if(file.exists(trackFile)) {
			# Precomputed depth track
			trk <- readRDT(trackFile)
		} else {
			# Run samtools
			command <- sprintf(
				"samtools depth -aa -q %i -Q %i -r chr%s:%i-%i \"%s\" | RLE",
				qBase, qMap, chrom, start, end, bamFile
			)
			tab <- read.table(
				pipe(command),
				sep="\t", header=FALSE, quote=NULL, comment.char="",
				col.names = c("chrom", "start", "end", "value"),
				colClasses = c("character", "integer", "integer", "integer")
			)
			
			# Depth track
			tab$name <- ""
			tab$chrom <- factor(sub("^chr", "", tab$chrom), chromosomes)
			tab$strand <- factor(NA, c("-","+"))
			trk <- track.table(
				tab,
				.name = sample,
				.organism = "Human",
				.assembly = assembly
			)
			
			# Export
			saveRDT(trk, file=trackFile)
		}
		
		return(trk)
	}

	# BAM file
	bamFile <- sprintf("%s/%s.DNA.MD.sort.bam", bamDir, sample)
	if(!file.exists(bamFile)) stop("\"", bamFile, "\" doesn't exist")
	
	# Annotation of the transcript of interest
	gene <- exons
	
	# Chromosome of the gene
	chrom <- unique(gene$chrom)
	if(length(chrom) != 1L) stop("Gene ", symbol, " on ", length(chrom), " chromosomes")
	
	# Exons (without redundancy)
	ano <- unique(gene[, c("start","end") ])
	
	# Genomic intervals (1-based, start & end included)
	breaks <- sort(unique(c(ano$start - 0.5, ano$end + 0.5)))
	ano <- data.frame(
		start = ceiling(head(breaks, -1)),
		end = floor(tail(breaks, -1))
	)
	rownames(ano) <- paste(ano$start, ano$end, sep="-")
	
	# Identify exons
	transcripts <- unique(gene$transcript)
	for(i in 1:length(transcripts)) {
		tab <- gene[ gene$transcript == transcripts[i] ,]
		for(j in 1:nrow(tab)) {
			# Genomic interval matching this exon
			match <- ano$start < tab[j,"end"] & ano$end > tab[j,"start"]
			ano[ match , transcripts[i] ] <- tab[j,"groupPosition"]
		}
	}
	
	# Sequencing depth in exons
	x <- which(apply(!is.na(ano[,-c(1,2),drop=FALSE]), 1, any))
	for(i in x) {
		trk <- depth(sample, bamFile, chrom=gene[1,"chrom"], start=ano[i,"start"], end=ano[i,"end"], trackDir=trackDir)
		ano[i,"depth"] <- sum(with(trk$extract(), (end-start+1L)*value)) / (ano[i,"end"] - ano[i,"start"] + 1L)
	}
	
	if(nrow(evt) > 0L) {
		# Normalized coordinates
		for(i in 1:nrow(evt)) {
			for(site in c("left", "right")) {
				# Overlapping feature
				if(chrom == evt[i, sprintf("%s.chrom", site) ]) {
					ovl <- which(ano$start <= evt[i, sprintf("%s.pos", site) ] & ano$end >= evt[i, sprintf("%s.pos", site) ])
					if(length(ovl) == 1L) {
						# Relative to the overlapped feature
						rel.pos <- as.integer(evt[i, sprintf("%s.pos", site) ] - ano[ovl,"start"])
						rel.size <- as.integer(ano[ovl,"end"] - ano[ovl,"start"])
						if(rel.size == 0L) { nrm <- ovl - 1L + 0.5
						} else             { nrm <- ovl - 1L + rel.pos / rel.size
						}
					} else if(length(ovl) > 1L) {
						stop("Ambiguity during feature overlap")
					} else if(all(evt[i, sprintf("%s.pos", site) ] < ano$start)) {
						# Before the gene
						nrm <- -0.5
					} else if(all(evt[i, sprintf("%s.pos", site) ] > ano$end)) {
						# After the gene
						nrm <- nrow(ano) + 0.5
					} else {
						stop("Unexpected case")
					}
				} else {
					# On another chromosome (fusion)
					nrm <- NA
				}
				
				# Store normalize coordinate
				evt[ i , sprintf("%s.nrm", site) ] <- nrm
			}
		}
		
		# Class color
		evt$color <- "red"
		evt$color[ evt$class == "annotated" ] <- "royalblue"
		evt$color[ evt$class == "plausible" ] <- "forestgreen"
		evt$color[ evt$class == "anchored" ] <- "orange"
		
		# Line
		evt$lty <- ifelse(evt$filter, "solid", "dotted")
		evt$lwd <- ifelse(evt$filter, 2, 1)
	}
	
	# Image file
	width <- 200 + nrow(ano) * 30
	height <- 460 + length(transcripts) * 40
	if(any(!is.na(evt$ID))) { IDs <- sprintf(" - %s", paste(na.omit(evt$ID[1:10]), collapse=" - "))
	} else                  { IDs <- ""
	}
	file <- sprintf("%s/%s - %s%s.png", outDir, symbol, sample, IDs)
	png(file=file, width=width, height=height, res=100)
	
	# Layout
	layout(matrix(1:3, ncol=1), heights=c(lcm(4), 1, lcm(4)))
	par(oma=c(3,1,3,1), cex=1)
	xlim <- c(-0.5, nrow(ano)+0.5)
	
	# Annotated junctions
	if(any(evt$class == "annotated")) { ymax <- log(max(evt[ evt$class == "annotated" , "reads" ]), 10)
	} else                            { ymax <- 1
	}
	par(mar=c(0,7,0,0))
	plot(x=NA, y=NA, xlim=xlim, ylim=c(0, ymax*1.05), xlab="", ylab="Reads", xaxs="i", xaxt="n", yaxt="n", yaxs="i", bty="n", las=2)
	for(i in which(evt$class == "annotated")) {
		# Coordinates
		x0 <- evt[i,"left.nrm"]
		x1 <- evt[i,"right.nrm"]
		y0 <- 0
		y1 <- log(evt[i,"reads"], 10)
		
		# Plot junction
		graphics::xspline(x=c(x0, x0, (x0+x1)/2, x1, x1), y=c(y0, y1, y1, y1, y0), shape=shape, lwd=evt[i,"lwd"], border=evt[i,"color"], col=evt[i,"color"], lty=evt[i,"lty"])
		
		# ID of candidates junctions
		if(!is.na(evt[i,"ID"])) text(x=(x0+x1)/2, y=y1, labels=evt[i,"ID"], col=evt[i,"color"], adj=c(0.5, -0.2), xpd=NA)
	}
	at <- 0:ceiling(ymax)
	axis(side=2, at=at, labels=10^at, las=2)
	at <- rep(2:9, ceiling(ymax)) * rep(10^(0:(ceiling(ymax)-1)), each=8)
	axis(side=2, at=log(at, 10), labels=FALSE, tcl=-0.3)
	
	# Main title
	title(main=sample, adj=0, outer=TRUE)
	title(main=symbol, adj=1, outer=TRUE)
	
	# Transcripts
	par(mar=c(1,7,1,0))
	plot(x=NA, y=NA, xlim=xlim, ylim=c(0, length(transcripts)), xlab="", ylab="", xaxs="i", xaxt="n", yaxt="n", yaxs="i", bty="n")
	for(i in 1:length(transcripts)) {
		# Included exons
		x <- ano[[ transcripts[i] ]]
		
		# Plot
		level <- ano[ which(!is.na(x)) , "depth" ] / max(ano$depth, na.rm=TRUE)
		level[ is.nan(level) ] <- 0L
		segments(x0=head(which(!is.na(x)), 1)-0.5, x1=tail(which(!is.na(x)), 1)-0.5, y0=i-0.5, y1=i-0.5)
		rect(xleft=which(!is.na(x))-1, xright=which(!is.na(x)), ybottom=i-0.9, ytop=i-0.1, col=grey(1 - level), border="#000000")
		text(x=which(!is.na(x))-0.5, y=i-0.5, adj=c(0.5, 0.5), labels=x[!is.na(x)], col=ifelse(level < 0.5, "black", "white"))
		mtext(side=2, at=i-0.5, text=sub(" .+$", "", transcripts[i]), las=2, line=0)
	}
	
	# Other junctions
	if(any(evt$class != "annotated")) { yaxt="s"; ymax <- log(max(evt[ evt$class != "annotated" , "reads" ]), 10)
	} else                            { yaxt="n"; ymax <- 1
	}
	par(mar=c(0,7,0,0))
	plot(x=NA, y=NA, xlim=xlim, ylim=c(ymax*1.05, 0), xlab="", ylab="Reads", xaxs="i", xaxt="n", yaxt="n", yaxs="i", bty="n", las=2)
	for(i in which(evt$class != "annotated")) {
		# Coordinates
		y0 <- 0
		y1 <- log(evt[i,"reads"], 10)
		
		left.genes <- strsplit(evt[i,"left.genes"], split=",")
		right.genes <- strsplit(evt[i,"right.genes"], split=",")
		is.left <- symbol %in% left.genes
		is.right <- symbol %in% right.genes
		
		if(is.left && !is.right) {
			# Fusion, right is partner
			x <- evt[i,"left.nrm"]
			if(!is.na(evt[i,"ID"])) {
				if(evt[i,"right.genes"] == "") { label <- evt[i,"ID"]
				} else                         { label <- sprintf("%s > %s", evt[i,"ID"], evt[i,"right.genes"])
				}
			}
			fusion <- TRUE
		} else if(!is.left && is.right) {
			# Fusion, left is partner
			x <- evt[i,"right.nrm"]
			if(!is.na(evt[i,"ID"])) {
				if(evt[i,"left.genes"] == "") { label <- evt[i,"ID"]
				} else                        { label <- sprintf("%s > %s", evt[i,"ID"], evt[i,"left.genes"])
				}
			}
			fusion <- TRUE
		} else {
			# Both are on the gene of interest
			fusion <- FALSE
		}
		
		if(fusion) {
			# Plot
			w <- (par("usr")[2] - par("usr")[1]) / 100
			graphics::segments(x0=x, x1=x, y0=y0, y1=y1, lwd=evt[i,"lwd"], col=evt[i,"color"], lty=evt[i,"lty"])
			graphics::segments(x0=x-0.3, x1=x+0.3, y0=y1, y1=y1, lwd=evt[i,"lwd"], col=evt[i,"color"], lty=evt[i,"lty"])
			
			# ID of candidates junctions
			if(!is.na(evt[i,"ID"])) text(x=x, y=y1, labels=label, col=evt[i,"color"], adj=c(0.5, 1.4), xpd=NA)
		} else {
			# Coordinates
			x0 <- evt[i,"left.nrm"]
			x1 <- evt[i,"right.nrm"]
			if(is.na(x0)) stop("Unexpectedly NA normalised coordinate (left)")
			if(is.na(x1)) stop("Unexpectedly NA normalised coordinate (right)")
			
			# Splicing event - plot
			graphics::xspline(x=c(x0, x0, (x0+x1)/2, x1, x1), y=c(y0, y1, y1, y1, y0), shape=shape, lwd=evt[i,"lwd"], border=evt[i,"color"], col=evt[i,"color"], lty=evt[i,"lty"])
			
			# ID of candidates junctions
			if(!is.na(evt[i,"ID"])) text(x=(x0+x1)/2, y=y1, labels=evt[i,"ID"], col=evt[i,"color"], adj=c(0.5, 1.2), xpd=NA)
		}
	}
	at <- 0:ceiling(ymax)
	axis(side=2, at=at, labels=10^at, las=2)
	at <- rep(2:9, ceiling(ymax)) * rep(10^(0:(ceiling(ymax)-1)), each=8)
	axis(side=2, at=log(at, 10), labels=FALSE, tcl=-0.3)
	
	# Legend
	legend(
		x="bottom", horiz=TRUE, inset=-0.3, xpd=NA, bty="n",
		lwd = c(2, 2, 2, 2, 1),
		lty = c("solid", "solid", "solid", "solid", "dotted"),
		col = c("royalblue", "forestgreen", "orange", "red", "black"),
		legend = c("annotated", "plausible", "anchored", "unknown", "filtered out")
	)
	
	void <- dev.off()
	
	invisible(TRUE)
	
}

# Converts an integer vector to another base, provided the digits to use
rebase <- function(x, base) {
	out <- character(length(x))
	if(length(x) > 0L) for(i in 1:length(x)) {
		n <- x[i]
		tmp <- integer(0)
		while(n > 0L) {
			tmp <- c(n %% length(base), tmp)
			n <- n %/% length(base)
		}
		if(length(tmp) == 0L) tmp <- 0L
		out[i] <- paste(base[tmp+1L], collapse="")
	}
	return(out)
}

# Export simplified table
exportCandidates <- function(events, groups, sites, I, S, events.filter.all, fusions, file="out/Candidates.tsv") {
	# Output column names
	columns <- c(
		"ID", "junction", "class", "recurrence", "sample", "reads", "fusion",
		"left.chrom", "left.pos", "left.genes", "left.exons.all", "left.transcripts.preferred", "left.exons.preferred", "left.depth", "left.PSI",
		"right.chrom", "right.pos",  "right.genes", "right.exons.all", "right.transcripts.preferred", "right.exons.preferred", "right.depth", "right.PSI"
	)
	
	# Events with at least 1 positive sample
	EOI <- rownames(events.filter.all)[ apply(events.filter.all, 1, any) ]
	if(length(EOI) > 0) {
		tab <- events[ EOI ,]
		
		# Add info on sites
		left <- sites[ sprintf("%s:%i", tab$left.chrom, tab$left.pos) ,]
		colnames(left) <- sprintf("left.%s", colnames(left))
		right <- sites[ sprintf("%s:%i", tab$right.chrom, tab$right.pos) ,]
		colnames(right) <- sprintf("right.%s", colnames(right))
		tab <- cbind(tab, left, right)
		
		# Gene-fusion or splicing event
		left.genes <- strsplit(tab$left.genes, split=",")
		right.genes <- strsplit(tab$right.genes, split=",")
		
		# Gene fusion or not
		fusion <- rep(TRUE, nrow(tab))
		fusion[ sapply(mapply(left.genes, right.genes, FUN=intersect), length) >= 1L ] <- FALSE
		fusion[ sapply(left.genes, length) == 0L & sapply(right.genes, length) == 0L ] <- FALSE
		fusion[ tab$left.chrom != tab$right.chrom ] <- TRUE
		fusion[ tab$right.pos - tab$left.pos > 500e3 ] <- TRUE
		tab$fusion <- fusion
		
		# Sequencing depth per sample
		depth <- I + S
		
		# Percentage Splice In
		PSI <- I / (I + S)
		
		# For each EOI, the list of samples passing filters (tab is pre-filtered so at least one)
		samples <- apply(events.filter.all[ rownames(tab) ,, drop=FALSE ], 1, function(x) { names(which(x)) })
		names(samples) <- NULL
		samples.n <- sapply(samples, length)
		samples <- unlist(samples)
		samples.i <- match(samples, colnames(I))
		
		# Storage (as a list, then a data.frame)
		out <- list()
		
		# EOI name, replicated for each sample
		out$junction <- rep(rownames(tab), samples.n)
		
		# Event row, replicated for each sample
		i <- rep(1:nrow(tab), samples.n)
		for(k in colnames(tab)) out[[k]] <- tab[i,k]
		
		# One row per sample
		out$sample <- samples
		
		# I, S (and derivatives) rows corresponding to left or right sites of each EOI
		groups$index <- 1:nrow(groups)
		groups.left.index <- with(groups[ groups$side == "left" ,], tapply(X=index, INDEX=event, FUN=unique))
		groups.right.index <- with(groups[ groups$side == "right" ,], tapply(X=index, INDEX=event, FUN=unique))
		
		# Replicate for each sample
		EOI.left  <- rep(groups.left.index[ rownames(tab) ], samples.n)
		EOI.right <- rep(groups.right.index[ rownames(tab) ], samples.n)
		
		# Event supporting reads, for each sample
		out$reads <- I[ cbind(EOI.left, samples.i) ]
		
		# PSI on each side, for each sample
		out$left.PSI  <- PSI[ cbind(EOI.left, samples.i) ]
		out$right.PSI <- PSI[ cbind(EOI.right, samples.i) ]
		
		# Depth on each side, for each sample
		out$left.depth  <- depth[ cbind(EOI.left, samples.i) ]
		out$right.depth <- depth[ cbind(EOI.right, samples.i) ]
		
		# Recurrence, replicated for each sample
		out$recurrence <- rep(samples.n, samples.n)
		
		# To data.frame
		if(length(unique(sapply(out, length))) != 1L) stop("Candidate gathering inconsistency")
		out <- as.data.frame(out)
		
		# Prioritize
		out <- out[ order(out$reads, decreasing=TRUE) ,]
		
		# Event ID
		event.ID <- factor(out$junction, levels=unique(out$junction))
		levels(event.ID) <- rebase((1:length(levels(event.ID)))-1L, LETTERS)
		
		# Sample ID
		regex <- "^.+_S([0-9])+$"
		samples <- unique(out$sample)
		if(all(grepl(regex, samples)) && !any(duplicated(as.integer(sub(regex, "\\1", samples))))) {
			# Illumina sample pattern : use sample sheet order
			sample.ID <- as.integer(sub(regex, "\\1", out$sample))
		} else {
			# No pattern : use alphabetical order
			sample.ID <- as.integer(factor(out$sample))
		}
		
		# Unique candidate ID
		ID <- paste(event.ID, sample.ID, sep="")
		if(any(duplicated(ID))) stop("Unicity error")
		out$ID <- ID
		
		# Filter and sort columns
		out <- out[, columns ]
		rownames(out) <- NULL
	} else {
		# Empty table
		out <- matrix(NA, nrow=0, ncol=length(columns), dimnames=list(NULL, columns))
	}
	
	# Filter out fusions
	if(!isTRUE(fusions)) out <- out[ !out$fusion ,]
	
	# Export
	if(!is.na(file)) write.table(out, file=file, row.names=FALSE, na="", sep="\t")
	
	invisible(out)
}



timedMessage("Loading dependencies...")

library(Rgb)
library(parallel)

timedMessage("Parsing symbol list...")

# Symbols of genes of interest
if(identical(symbols, "all")) {
	# All symbols defined in the genome
	symbolList <- NULL
} else if(identical(symbols, "target")) {
	# All symbols defined in the target GTF
	gtf <- read.gtf(pipe(sprintf("awk '$3 == \"gene\" { print }' \"%s\"", targetFile)))
	if("gene" %in% colnames(gtf))           { symbolList <- gtf$gene
	} else if("gene_id" %in% colnames(gtf)) { symbolList <- gtf$gene_id
	} else                                  { stop("Was expecting a 'gene' or 'gene_id' feature in targetGTF")
	}
} else {
	# User-provided list of symbols
	symbolList <- strsplit(symbols, split=",")[[1]]
}

timedMessage("Parsing annotation...")

exons <- readRDT(exonFile)

timedMessage("Parsing collected events...")

I <- readRDS("I.rds")
S <- readRDS("S.rds")
groups <- readRDS("groups.rds")
sites <- readRDS("sites.rds")
events <- readRDS("events.rds")

timedMessage("Filtering classes...")

events <- filterClass(events, classes=classList)

timedMessage("Filtering symbols...")

out <- filterSymbol(events, groups, sites, symbols=symbolList)
events <- out$events
groups <- out$groups
sites <- out$sites

timedMessage("Filtering I...")

events.filter.I <- filterI(events, I, min.I)

timedMessage("Filtering PSI...")

events.filter.PSI <- filterPSI(events, groups, I, S, min.PSI)

timedMessage("Overlapping filters...")

# 'event' x 'sample' passing all filters at once
events.filter.all <- events$filter.class & events$filter.symbol & events.filter.I & events.filter.PSI

# Stats
message("- ", sum(events.filter.all), " positives in ", sum(apply(events.filter.all, 1, any)), " junctions of interest")

timedMessage("Exporting candidates...")

# Output directory
outDir <- sprintf(
	"I-%i_PSI-%g_%s(%i)_%s_%s_%s",
	min.I,
	min.PSI,
	substr(symbols, 1, 50),
	length(strsplit(symbols, split=",")[[1]]),
	classes,
	gsub(":", "-", focus),
	ifelse(fusions, "fusions", "no-fusions")
)
dir.create(outDir)

# Candidates
candidates <- exportCandidates(events, groups, sites, I, S, events.filter.all, fusions, file=sprintf("%s/Candidates.tsv", outDir))

# Sequencing depth data directory
dir.create("depth")

if(isTRUE(plot) && nrow(candidates) > 0L) {
	# Output directory
	dir.create(sprintf("%s/plots", outDir))
	
	# Pre-filter data to minimize transfers during parallelization
	timedMessage("Pre-filtering data for plots...")
	toPlot <- preparePlots(candidates, events, groups, sites, events.filter.all)
	
	# Create a cluster for parallelization
	timedMessage("Plotting on ", ncores, " CPUs...")
	cluster <- makeCluster(spec=ncores)
	void <- clusterMap(
		cl = cluster,
		fun = plot.normalized,
		evt = toPlot$events,
		sample = toPlot$sample,
		symbol = toPlot$symbol,
		exons = toPlot$exons,
		MoreArgs = list(
			outDir = sprintf("%s/plots", outDir),
			bamDir = ".",
			trackDir = "depth"
		)
	)
	stopCluster(cluster)
	
} else timedMessage("Nothing to plot")

timedMessage("done")

