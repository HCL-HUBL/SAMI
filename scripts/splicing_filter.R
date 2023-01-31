#!/usr/bin/env Rscript

# Collect CLI arguments
args <- commandArgs(TRUE)
if(length(args) != 8L) stop("USAGE : ./splicing_filter.R NCORES exons.rdt PLOT MIN_I MIN_PSI SYMBOLS|all CLASSES FOCUS")
ncores <- as.integer(args[1])
exonFile <- args[2]
plot <- as.logical(args[3])
min.I <- as.integer(args[4])
min.PSI <- as.double(args[5])
symbols <- args[6]
if(identical(symbols, "all")) { symbolList <- NULL
} else                        { symbolList <- strsplit(symbols, split=",")[[1]]
}
if(length(symbolList) > 10L) symbols <- sprintf("%i-symbols", length(symbolList))
classes <- args[7]
classList <- strsplit(classes, split=",")[[1]]
focus <- args[8]
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

# Filter events with both sites overlapping a gene of interest
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
	
	# Event level (left AND right)
	events.filter.symbol <- tapply(X=groups$filter.symbol, INDEX=groups$event, FUN=all)
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
	left  <- groups.filter.PSI[ groups$side == "left" ,][ rownames(events) ,]
	right <- groups.filter.PSI[ groups$side == "right" ,][ rownames(events) ,]
	events.filter.PSI <- left & right
	
	return(events.filter.PSI)
}

# Pre-filter data to minimize transfers during parallelization
preparePlots <- function(candidates, events, groups, sites, events.filter.all) {
	# All genes involved for each candidate
	left.genes <- strsplit(candidates$left.genes, ",")
	right.genes <- strsplit(candidates$right.genes, ",")
	genes <- mapply(left.genes, right.genes, FUN=unique)
	
	# All plots to produce
	n <- sapply(genes, length)
	toPlot <- data.frame(
		### junction = rep(candidates$junction, n),
		symbol = sub(" \\([+-]\\)$", "", unlist(genes)),
		sample = rep(candidates$sample, n)
	)
	
	# Reduce redundancy
	toPlot <- unique(toPlot)
	toPlot <- as.list(toPlot)
	
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
				ovl <- which(ano$start <= evt[i,site] & ano$end >= evt[i,site])
				if(length(ovl) == 1L) {
					# Relative to the overlapped feature
					nrm <- ovl - 1L + (evt[i,site] - ano[ovl,"start"]) / (ano[ovl,"end"] - ano[ovl,"start"])
				} else if(length(ovl) > 1L) {
					stop("Ambiguity during feature overlap")
				} else if(all(evt[i,site] < ano$start)) {
					# Before the gene
					nrm <- -0.5
				} else if(all(evt[i,site] > ano$end)) {
					# After the gene
					nrm <- nrow(ano) + 0.5
				} else {
					stop("Unexpected case")
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
	plot(x=NA, y=NA, xlim=xlim, ylim=c(0, ymax), xlab="", ylab="Reads", xaxs="i", xaxt="n", yaxt="n", yaxs="i", bty="n", las=2)
	for(i in which(evt$class == "annotated")) {
		# Coordinates
		x0 <- evt[i,"left.nrm"]
		x1 <- evt[i,"right.nrm"]
		y0 <- 0
		y1 <- log(evt[i,"reads"], 10)
		
		# Plot junction
		graphics::xspline(x=c(x0, x0, (x0+x1)/2, x1, x1), y=c(y0, y1, y1, y1, y0), shape=shape, lwd=evt[i,"lwd"], border=evt[i,"color"], col=evt[i,"color"], lty=evt[i,"lty"], xpd=NA)
		
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
	plot(x=NA, y=NA, xlim=xlim, ylim=c(ymax, 0), xlab="", ylab="Reads", xaxs="i", xaxt="n", yaxt="n", yaxs="i", bty="n", las=2)
	for(i in which(evt$class != "annotated")) {
		# Coordinates
		x0 <- evt[i,"left.nrm"]
		x1 <- evt[i,"right.nrm"]
		y0 <- 0
		y1 <- log(evt[i,"reads"], 10)
		
		# Plot junction
		graphics::xspline(x=c(x0, x0, (x0+x1)/2, x1, x1), y=c(y0, y1, y1, y1, y0), shape=shape, lwd=evt[i,"lwd"], border=evt[i,"color"], col=evt[i,"color"], lty=evt[i,"lty"], xpd=NA)
		
		# ID of candidates junctions
		if(!is.na(evt[i,"ID"])) text(x=(x0+x1)/2, y=y1, labels=evt[i,"ID"], col=evt[i,"color"], adj=c(0.5, 1.2), xpd=NA)
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
exportCandidates <- function(events, groups, sites, I, S, events.filter.all, file="out/Candidates.csv") {
	# Output column names
	columns <- c(
		"ID", "junction", "chrom", "class", "recurrence", "sample", "reads",
		"left", "left.genes", "left.exons.all", "left.transcripts.preferred", "left.exons.preferred", "left.depth", "left.PSI",
		"right", "right.genes", "right.exons.all", "right.transcripts.preferred", "right.exons.preferred", "right.depth", "right.PSI"
	)
	
	# Events with at least 1 positive sample
	EOI <- rownames(events.filter.all)[ apply(events.filter.all, 1, any) ]
	if(length(EOI) > 0) {
		tab <- events[ EOI ,]
		
		# Add info on sites
		left <- sites[ sprintf("%s:%i", tab$chrom, tab$left) ,]
		colnames(left) <- sprintf("left.%s", colnames(left))
		right <- sites[ sprintf("%s:%i", tab$chrom, tab$right) ,]
		colnames(right) <- sprintf("right.%s", colnames(right))
		tab <- cbind(tab, left, right)
		
		# Sequencing depth per sample
		depth <- I + S
		
		# Percentage Splice In
		PSI <- I / (I + S)
		
		# Duplicate per sample
		out <- vector(mode="list", nrow(tab))
		for(i in 1:nrow(tab)) {
			
			if(i %% 500L == 0L) message("- ", i, "/", nrow(tab))
		
			# Positive samples for this event
			samples <- names(which(events.filter.all[ rownames(tab)[i] ,]))
			
			# Add sample-level data
			is.event <- groups$event == rownames(tab)[i]
			is.left  <- is.event & groups$side == "left"
			is.right <- is.event & groups$side == "right"
			out[[i]] <- data.frame(
				junction = rownames(tab)[i],
				tab[i,],
				sample = samples,
				reads = I[ is.left , samples ],
				left.PSI  = PSI[ is.left , samples ],
				right.PSI = PSI[ is.right , samples ],
				left.depth  = depth[ is.left , samples ],
				right.depth = depth[ is.right , samples ],
				recurrence = length(samples),
				row.names = NULL
			)
		}
		
		# Merge events
		tab <- do.call(rbind, out)
		
		# Prioritize
		tab <- tab[ order(tab$reads, decreasing=TRUE) ,]
		
		# Event ID
		event.ID <- factor(tab$junction, levels=unique(tab$junction))
		levels(event.ID) <- rebase((1:length(levels(event.ID)))-1L, LETTERS)
		
		# Sample ID
		if(all(grepl("^.+_S([0-9])+$", tab$sample))) {
			# Illumina sample pattern : use sample sheet order
			sample.ID <- as.integer(sub("^.+_S([0-9]+)$", "\\1", tab$sample))
		} else {
			# No pattern : use alphabetical order
			sample.ID <- as.integer(factor(tab$sample))
		}
		
		# Unique candidate ID
		ID <- paste(event.ID, sample.ID, sep="")
		if(any(duplicated(ID))) stop("Unicity error")
		tab$ID <- ID
		
		# Filter and sort columns
		tab <- tab[, columns ]
		rownames(tab) <- NULL
	} else {
		# Empty table
		tab <- matrix(NA, nrow=0, ncol=length(columns), dimnames=list(NULL, columns))
	}
	
	# Export
	if(!is.na(file)) write.csv2(tab, file=file, row.names=FALSE, na="")
	
	invisible(tab)
}



timedMessage("Loading dependencies...")

library(Rgb)
library(openxlsx)
library(parallel)

### timedMessage("Parsing preferred transcript file...")
### preferred <- parsePreferred(transcriptFile)

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
outDir <- sprintf("I-%i_PSI-%g_%s_%s_%s", min.I, min.PSI, symbols, classes, gsub(":", "-", focus))
dir.create(outDir)

# Candidates
candidates <- exportCandidates(events, groups, sites, I, S, events.filter.all, file=sprintf("%s/Candidates.csv", outDir))

if(isTRUE(plot) && nrow(candidates) > 0L) {
	
	# Output directory
	dir.create("depth")
	dir.create(sprintf("%s/plots", outDir))
	
	# Pre-filter data to minimize transfers during parallelization
	timedMessage("Pre-filtering data for plots...")
	toPlot <- preparePlots(candidates, events, groups, sites, events.filter.all)
	
	# Create a cluster for parallelization
	timedMessage("Plotting ", length(toPlot$sample), " genes & samples on ", ncores, " CPUs...")
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

