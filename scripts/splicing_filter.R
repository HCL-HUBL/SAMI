#!/usr/bin/env Rscript

# Collect CLI arguments
### args <- commandArgs(TRUE)
args <- c(10, "../../../store/exons.GRCh38.rdt", FALSE, FALSE, 30, 0.1, "all", "plausible", "none")
if(length(args) != 9L) stop("USAGE : ./splicing_filter.R NCORES exons.rdt XLSX PLOT MIN_I MIN_PSI SYMBOLS|all CLASSES FOCUS")
ncores <- as.integer(args[1])
exonFile <- args[2]
xlsx <- as.logical(args[3])
plot <- as.logical(args[4])
min.I <- as.integer(args[5])
min.PSI <- as.double(args[6])
symbols <- args[7]
if(identical(symbols, "all")) { symbolList <- NULL
} else                        { symbolList <- strsplit(symbols, split=",")[[1]]
}
if(length(symbolList) > 10L) symbols <- sprintf("%i-symbols", length(symbolList))
classes <- args[8]
classList <- strsplit(classes, split=",")[[1]]
focus <- args[9]
if(identical(focus, "none")) { focusList <- NULL
} else                       { focusList <- strsplit(focus, split=",")[[1]]
}



# Print a log message with date and time
timedMessage <- function(...) {
	message(Sys.time(), " : ", ...)
}

# Similar to lapply() but apparently faster
loop <- function(X, FUN, ...) {
	out <- list()
	for(i in 1:length(X)) {
		out[[i]] <- FUN(X[[i]], ...)
	}
	return(out)
}

# Aggregate for each target (candidate event of interest) junctions sharing the same left or right boundary (results in event duplication)
extendEvents <- function(events, groups) {
	# Event list for each site
	sites.events <- tapply(X=groups$event, INDEX=groups$site, FUN=c)
	
	# IDs of the two sites of each target event
	targets.left  <- sprintf("%s:%i", events$chrom, events$left)
	targets.right <- sprintf("%s:%i", events$chrom, events$right)
	
	# Aggregate events sharing the same left as the target event
	targets.lefts <- sites.events[ targets.left ]
	names(targets.lefts) <- rownames(events)
	targets.lefts <- data.frame(
		target = rep(names(targets.lefts), sapply(targets.lefts, length)),
		event = unlist(targets.lefts),
		side = "left"
	)
	
	# Aggregate events sharing the same right as the target event
	targets.rights <- sites.events[ targets.right ]
	names(targets.rights) <- rownames(events)
	targets.rights <- data.frame(
		target = rep(names(targets.rights), sapply(targets.rights, length)),
		event = unlist(targets.rights),
		side = "right"
	)
	
	# Collect left and right events
	targets <- rbind(targets.lefts, targets.rights)
	targets <- targets[ order(targets$target, targets$event) ,]
	rownames(targets) <- NULL
	
	return(targets)
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
	
	# Filter status
	evt$filter <- evt[, sprintf("filter.%s", sample) ]
	
	# Supporting reads
	evt$reads <- evt[, sprintf("I.%s", sample) ]
	e <- evt[ evt$reads > 0L ,]
	
	if(nrow(e) > 0L) {
		# Normalized coordinates
		for(i in 1:nrow(e)) {
			for(site in c("left", "right")) {
				# Overlapping feature
				ovl <- which(ano$start <= e[i,site] & ano$end >= e[i,site])
				if(length(ovl) == 1L) {
					# Relative to the overlapped feature
					nrm <- ovl - 1L + (e[i,site] - ano[ovl,"start"]) / (ano[ovl,"end"] - ano[ovl,"start"])
				} else if(length(ovl) > 1L) {
					stop("Ambiguity during feature overlap")
				} else if(all(e[i,site] < ano$start)) {
					# Before the gene
					nrm <- -0.5
				} else if(all(e[i,site] > ano$end)) {
					# After the gene
					nrm <- nrow(ano) + 0.5
				} else {
					stop("Unexpected case")
				}
				
				# Store normalize coordinate
				e[ i , sprintf("%s.nrm", site) ] <- nrm
			}
		}
		
		# Class color
		e$color <- "red"
		e$color[ e$class == "annotated" ] <- "royalblue"
		e$color[ e$class == "plausible" ] <- "forestgreen"
		e$color[ e$class == "anchored" ] <- "orange"
		
		# Line
		e$lty <- ifelse(e$filter, "solid", "dotted")
		e$lwd <- ifelse(e$filter, 2, 1)
	}
	
	# Image file
	width <- 200 + nrow(ano) * 30
	height <- 460 + length(transcripts) * 40
	file <- sprintf("%s/%s - %s.png", outDir, symbol, sample)
	png(file=file, width=width, height=height, res=100)

	# Layout
	layout(matrix(1:3, ncol=1), heights=c(lcm(4), 1, lcm(4)))
	par(oma=c(3,1,3,1), cex=1)
	xlim <- c(-0.5, nrow(ano)+0.5)
	
	# Annotated junctions
	if(any(e$class == "annotated")) { ymax <- log(max(e[ e$class == "annotated" , "reads" ]), 10)
	} else                          { ymax <- 1
	}
	par(mar=c(0,7,0,0))
	plot(x=NA, y=NA, xlim=xlim, ylim=c(0, ymax), xlab="", ylab="Reads", xaxs="i", xaxt="n", yaxt="n", yaxs="i", bty="n", las=2)
	for(i in which(e$class == "annotated")) {
		# Coordinates
		x0 <- e[i,"left.nrm"]
		x1 <- e[i,"right.nrm"]
		y0 <- 0
		y1 <- log(e[i,"reads"], 10)
		
		# Plot junction
		graphics::xspline(x=c(x0, x0, (x0+x1)/2, x1, x1), y=c(y0, y1, y1, y1, y0), shape=shape, lwd=e[i,"lwd"], border=e[i,"color"], col=e[i,"color"], lty=e[i,"lty"], xpd=NA)
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
	if(any(e$class != "annotated")) { yaxt="s"; ymax <- log(max(e[ e$class != "annotated" , "reads" ]), 10)
	} else                          { yaxt="n"; ymax <- 1
	}
	par(mar=c(0,7,0,0))
	plot(x=NA, y=NA, xlim=xlim, ylim=c(ymax, 0), xlab="", ylab="Reads", xaxs="i", xaxt="n", yaxt="n", yaxs="i", bty="n", las=2)
	for(i in which(e$class != "annotated")) {
		# Coordinates
		x0 <- e[i,"left.nrm"]
		x1 <- e[i,"right.nrm"]
		y0 <- 0
		y1 <- log(e[i,"reads"], 10)
		
		# Plot junction
		graphics::xspline(x=c(x0, x0, (x0+x1)/2, x1, x1), y=c(y0, y1, y1, y1, y0), shape=shape, lwd=e[i,"lwd"], border=e[i,"color"], col=e[i,"color"], lty=e[i,"lty"], xpd=NA)
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
		legend = c("annotated", "plausible", "anchored", "unknown", "filtered")
	)
	
	void <- dev.off()
	
	invisible(TRUE)
	
}

# Export detailed table in CSV
exportDetails <- function(tab, file="out/Details.csv") {
	# Columns
	col.now <- c("target", "site", "label", "chrom", "left", "right", "ID.symbol", "class", "exon.site", "exon.partner")
	col.clean <- c("Jonction d'intérêt", "Site d'épissage", "Alternative", "Chrom", "Gauche", "Droite", "Gene", "Classe", "Exon (site)", "Exon (partenaire)")
	
	if(!is.data.frame(tab)) {
		# Empty table
		exp <- matrix(NA, nrow=0, ncol=length(col.clean), dimnames=list(NULL, col.clean))
	} else {
		# Column selection
		exp <- tab[, col.now ]
		colnames(exp) <- col.clean
		
		# Rename I and round PSI
		for(sample in sub("^PSI\\.", "", grep("^PSI", colnames(tab), value=TRUE))) {
			exp[[ sprintf("%s.I", sample) ]] <- tab[[ sprintf("I.%s", sample) ]]
			exp[[ sprintf("%s.PSI", sample) ]] <- round(tab[[ sprintf("PSI.%s", sample) ]], 3)
		}
	}
	
	# Export
	if(!is.na(file)) write.csv2(exp, file=file, row.names=FALSE, na="")
	
	invisible(exp)
}

# Export detailed table in XLSX
formatDetails <- function(details, tab, file="out/Details.xlsx") {
	# Create workbook
	wb <- createWorkbook()
	addWorksheet(wb, "Junctions")
	
	timedMessage("- Printing data blocks...")
	
	if(nrow(details) == 0L) {
		# Empty table
		writeData(wb=wb, sheet=1L, x=details, borders="surrounding")
	} else {
		# Print data
		startRow <- 1L
		for(target in unique(details$"Jonction d'intérêt")) {
			
			message("-- ", target)
			
			# Add one target at the time with borders
			tmp <- details[ details$"Jonction d'intérêt" == target ,]
			tmp[ is.na(tmp) | tmp == 0L ] <- ""
			writeData(wb=wb, sheet=1L, startCol=1L, startRow=startRow, colNames=(startRow == 1L), x=tmp, borders="surrounding")
			startRow <- startRow + nrow(tmp) + (startRow == 1L)
			
			# Merge full subset height
			if(nrow(tmp) > 1L) {
				range <- (startRow - nrow(tmp)) : (startRow - 1L)
				for(k in match(c("Jonction d'intérêt", "Chrom", "Gene"), colnames(details))) mergeCells(wb, sheet=1L, cols=k, rows=range)
			}
			
			# Merge half subset height
			for(site in unique(details$"Site d'épissage")) {
				i <- which(details$"Jonction d'intérêt" == target & details$"Site d'épissage" == site) + 1L
				if(length(i) > 1L) {
					range <- range(i)
					for(k in match(c("Site d'épissage", "Exon (site)"), colnames(details))) mergeCells(wb, sheet=1L, cols=k, rows=range)
				}
			}
		}
		
		timedMessage("- Adding 'selected' style...")
		
		# Cell styles
		styles <- list(
			"selected"  = createStyle(fgFill="lightgrey"),
			"annotated" = createStyle(fgFill="lightblue"),
			"plausible" = createStyle(fgFill="lightgreen"),
			"anchored"  = createStyle(fgFill="yellow"),
			"unknown"   = createStyle(fgFill="pink")
		)
		
		# Add style to alternative A
		addStyle(
			wb=wb, sheet=1L, style=styles[["selected"]], stack=TRUE,
			rows = which(details$"Alternative" == "A") + 1L,
			cols = c(3, 5, 6, 8, 10:ncol(details)),
			gridExpand = TRUE
		)
		
		# Add style to cells of interest
		significant <- tab[, grep("^filter\\.", colnames(tab)) ]
		highlight <- list()
		for(class in c("annotated", "plausible", "anchored", "unknown")) {
			
			timedMessage("- Adding '", class, "' style...")
			
			# Cells of interest
			highlight <- significant & details$Classe == class
			
			# Add style
			addStyle(
				wb=wb, sheet=1L, style=styles[[class]], stack=TRUE,
				rows = c(
					row(highlight)[ highlight ] + 1L,
					row(highlight)[ highlight ] + 1L
				),
				cols = c(
					col(highlight)[ highlight ] * 2 + 9,
					col(highlight)[ highlight ] * 2 + 10
				)
			)
		}
	}
	
	timedMessage("- Writting file...")
	
	# Save workbook
	saveWorkbook(wb, file, overwrite=TRUE)
}

# Export simplified table
exportCandidates <- function(tab, file="out/Candidates.csv") {
	# Sample list
	samples <- sub("^filter\\.", "", grep("^filter\\.", colnames(tab), value=TRUE))
	
	# Process samples one at a time
	cand <- list()
	for(sample in samples) {
		# Junctions of interest for this sample
		jun <- unique(tab[ tab[[ sprintf("filter.%s", sample) ]] & tab$label == "A" , "target" ])
		if(length(jun) > 0L) {
			# Cells of interest for this sample
			rows <- which(tab$target %in% jun & tab$label == "A")
			tmp <- tab[ rows , c("target", "site", "chrom", "left", "right", "ID.symbol", "class", "exons.site", "exons.partner", "exon.transcript", "exon.site", "exon.partner", sprintf("PSI.%s", sample), sprintf("I.%s", sample), sprintf("depth.%s", sample)) ]
			tmp$is.left <- tmp$site == sprintf("%s:%i", tmp$chrom, tmp$left)
			tmp$site <- NULL
			
			# Merge row pairs
			mrg <- tmp[ tmp$is.left ,]
			mrg$is.left <- NULL
			mrg$sample <- sample
			mrg$reads <- mrg[[ sprintf("I.%s", sample) ]]
			mrg$PSI.left <- tmp[  tmp$is.left , sprintf("PSI.%s", sample) ]
			mrg$PSI.right <- tmp[ !tmp$is.left , sprintf("PSI.%s", sample) ]
			colnames(mrg)[ colnames(mrg) == "exon.site" ] <- "exon.left"
			colnames(mrg)[ colnames(mrg) == "exon.partner" ] <- "exon.right"
			colnames(mrg)[ colnames(mrg) == "exons.site" ] <- "exons.left"
			colnames(mrg)[ colnames(mrg) == "exons.partner" ] <- "exons.right"
			mrg$depth.left <- tmp[  tmp$is.left , sprintf("depth.%s", sample) ]
			mrg$depth.right <- tmp[ !tmp$is.left , sprintf("depth.%s", sample) ]
			mrg[[ sprintf("PSI.%s", sample) ]] <- NULL
			mrg[[ sprintf("I.%s", sample) ]] <- NULL
			mrg[[ sprintf("depth.%s", sample) ]] <- NULL
			
			# Aggregate
			cand[[ sample ]] <- mrg
		}
	}
	
	# Merge samples
	cand <- do.call(rbind, cand)
	
	# Intra-run recurrence, post-filtering
	cand$recurrence <- as.integer(table(cand$target)[ cand$target ])
	
	# Columns
	col.now <- c("target", "chrom", "left", "right", "ID.symbol", "class", "exons.left", "exons.right", "exon.transcript", "exon.left", "exon.right", "sample", "reads", "PSI.left", "PSI.right", "depth.left", "depth.right", "recurrence")
	col.clean <- c("Jonction d'intérêt", "Chrom", "Gauche", "Droite", "Gene", "Classe", "Exons (gauche)", "Exons (droite)", "Transcrit préféré", "Exon (gauche)", "Exon (droite)", "Echantillon", "Reads", "PSI (gauche)", "PSI (droite)", "Profondeur (gauche)", "Profondeur (droite)", "Récurrence intra-run")
	
	if(length(cand) == 0) {
		# Empty table
		cand <- matrix(NA, nrow=0, ncol=length(col.clean), dimnames=list(NULL, col.clean))
	} else {
		# Prioritize
		cand <- cand[ order(cand$reads, decreasing=TRUE) ,]
		
		# Rename columns
		cand <- cand[, col.now ]
		colnames(cand) <- col.clean
	}
	
	# Export
	if(!is.na(file)) write.csv2(cand, file=file, row.names=FALSE, na="")
	
	invisible(cand)
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

events$filter.class <- events$class %in% classList

timedMessage("Filtering symbols...")

# Symbol filtering in 'sites' (computation)
sites$filter.symbol <- filterSites(sites, symbols.filter=symbolList)

# Symbol filtering in 'groups' (propagation)
mrg <- merge(x=groups, y=sites[,"filter.symbol",drop=FALSE], by.x="site", by.y="row.names", all=TRUE, sort=FALSE)
if(!identical(mrg[,1:3], groups[,1:3])) stop("Merging error")
groups <- mrg

# Symbol filtering in 'events' (left AND right)
events.filter.symbol <- tapply(X=groups$filter.symbol, INDEX=groups$event, FUN=all)
events[ , "filter.symbol" ] <- as.logical(events.filter.symbol[ rownames(events) ])

timedMessage("Filtering I...")

# Deduplicate I at 'event' x 'sample' level
events.I <- I[ rownames(events) ,]

# Events with at least min.I supporting reads for each sample
events.filter.I <- events.I >= min.I

timedMessage("Filtering PSI...")

# PSI at 'groups' x 'sample' level
PSI <- I / (I + S)

# PSI filtering at 'groups' x 'sample' level
groups.filter.PSI <- !is.na(PSI) & PSI >= min.PSI

# PSI filtering at 'event' x 'sample' level (left AND right)
left  <- groups.filter.PSI[ groups$side == "left" ,][ rownames(events) ,]
right <- groups.filter.PSI[ groups$side == "right" ,][ rownames(events) ,]
events.filter.PSI <- left & right

timedMessage("Overlapping filters...")

# 'event' x 'sample' passing all filters at once
events.filter.all <- events$filter.class & events$filter.symbol & events.filter.I & events.filter.PSI

# Stats
message("- ", sum(events.filter.all), " positives in ", sum(apply(events.filter.all, 1, any)), " junctions of interest")



stop("That's all folks")

timedMessage("Exporting candidates...")

# Work in progress : easy merging
exportCandidates(out, file=sprintf("%s/Candidates.csv", outDir))





timedMessage("Extending events...")

targets <- extendEvents(events, groups)







timedMessage("Preparing output directory...")

outDir <- sprintf("I-%i_PSI-%g_%s_%s_%s", min.I, min.PSI, symbols, classes, gsub(":", "-", focus))
dir.create(outDir)
dir.create("depth")

if(isTRUE(plot) && is.data.frame(toPlot) && nrow(toPlot) > 0L) {
	
	timedMessage("Preparing genes to plot...")
	
	# Output directory
	dir.create(sprintf("%s/plots", outDir))
	
	# Pre-filter exons to minimize transfers during parallelization
	toPlot.exons <- list()
	for(symbol in unique(toPlot$symbol)) {
		gene <- exons$extract(exons$extract(,"symbol") == symbol)
		for(i in which(toPlot$symbol == symbol)) toPlot.exons[[i]] <- gene
	}
	
	# Pre-filter events to minimize transfers during parallelization
	toPlot.events <- list()
	for(symbol in unique(toPlot$symbol)) {
		match.symbol <- sapply(
			strsplit(gsub(" \\([+-]\\)", "", events$ID.symbol), split=", "),
			`%in%`, x=symbol
		)
		evt <- events[ match.symbol ,]
		for(i in which(toPlot$symbol == symbol)) toPlot.events[[i]] <- evt
	}
	
	timedMessage("Plotting ", nrow(toPlot), " genes & samples on ", ncores, " CPUs...")
	void <- clusterMap(
		cl = cluster,
		fun = plot.normalized,
		sample = toPlot$sample,
		symbol = toPlot$symbol,
		exons = toPlot.exons,
		evt = toPlot.events,
		MoreArgs = list(
			outDir = sprintf("%s/plots", outDir),
			bamDir = ".",
			trackDir = "depth"
		)
	)
	
} else timedMessage("Nothing to plot")

# Close the cluster
stopCluster(cluster)

timedMessage("Exporting details (CSV)...")

details <- exportDetails(out, file=sprintf("%s/Details.csv", outDir))

if(xlsx) {
	timedMessage("Exporting details (XLSX)...")
	formatDetails(details, out, file=sprintf("%s/Details.xlsx", outDir))
}

timedMessage("done")
