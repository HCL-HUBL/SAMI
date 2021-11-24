#!/usr/bin/env Rscript

# Collect CLI arguments
args <- commandArgs(TRUE)
if(length(args) != 8L) stop("USAGE : ./splicing_filter.R events.rds exons.rdt PLOT MIN_I MIN_PSI SYMBOLS|all CLASSES FOCUS")
eventFile <- args[1]
exonFile <- args[2]
plot <- as.logical(args[3])
min.I <- as.integer(args[4])
min.PSI <- as.double(args[5])
symbols <- args[6]
if(identical(symbols, "all")) { symbolList <- NULL
} else                        { symbolList <- strsplit(symbols, split=",")[[1]]
}
classes <- args[7]
classList <- strsplit(classes, split=",")[[1]]
focus <- args[8]
if(identical(focus, "none")) { focusList <- NULL
} else                       { focusList <- strsplit(focus, split=",")[[1]]
}



# Determine for each junction event whether it passes evidence strength filters or not
filterEvents <- function(events, min.PSI=0.01, min.I=3L) {
	# Samples / junctions with sufficient data
	isSignificant <- events[, grep("^PSI\\.", colnames(events)) ] >= min.PSI & events[, grep("^I\\.", colnames(events)) ] >= min.I
	colnames(isSignificant) <- sub("^PSI\\.", "filter.", colnames(isSignificant))
	events <- cbind(events, isSignificant)
	
	# Any sample passes filter for this junction
	events$anySignificant.filter <- apply(isSignificant, 1, any)
	
	return(events)
}

# Aggregate for each site junctions sharing the same left or right boundary (results in event duplication)
extendEvents <- function(events) {
	# Split event table by site (named list)
	sites <- split(events, events$site)
	
	# All considered junctions
	jun <- unique(events$ID)
	jun.left  <- sub("^([0-9XY]+):([0-9]+).([0-9]+|nosplice)$", "\\1:\\2", jun)
	jun.right <- sub("^([0-9XY]+):([0-9]+).([0-9]+|nosplice)$", "\\1:\\3", jun)
	
	# Group events sharing the same left or right boundary than the considered site
	left  <- sites[ jun.left ]
	right <- sites[ jun.right ]
	
	# Merge left and right collections site by site
	junctions <- mapply(left, right, FUN=rbind, SIMPLIFY=FALSE)
	names(junctions) <- jun
	
	return(junctions)
}

# Identify groups of events (1 junction of interest + alternatives) passing filters
filterExtended <- function(extended, symbols.filter=NULL, A.types.filter=c("unknown", "anchored", "plausible"), focus=NULL) {
	if(length(focus) > 0L) {
		# Provided list of junctions to keep
		missing <- setdiff(focus, names(extended))
		if(length(missing) > 0L) stop("Requested junction(s) missing : ", paste(missing, collapse=", "))
		keep <- names(extended) %in% focus
		names(keep) <- names(extended)
	} else {
		# Junctions of interest
		keep <- sapply(
			extended,
			function(x) {
				# At least one of the two sites is of interest
				out <- TRUE
				
				# Genes of interest
				if(length(symbols.filter) > 0L) out <- out && any(x$ID.symbol %in% symbols.filter)
				
				# A is significant
				A <- duplicated(x$ID) | duplicated(x$ID, fromLast=TRUE)
				out <- out && all(x[ A , "anySignificant.filter" ])
				
				# A type is of interest
				out <- out && all(x[ A , "class" ] %in% A.types.filter)
				
				return(out)
			}
		)
	}
	
	return(keep)
}

# Determine whether the considered genomic position corresponds to an exon boundary or not
annotateSplicingSite <- function(chrom, position, exons) {
	anno <- NULL
	
	# Is on correct chromosome
	isOnChrom <- exons$extract(,"chrom") == chrom
	
	# Position of interest is exon start
	i <- which(isOnChrom & exons$extract(,"start") == position)
	if(length(i) > 0L) {
		anno <- rbind(
			anno,
			unique(with(exons$extract(i), data.frame(exon=groupPosition, end="left")))
		)
	}
	
	# Position of interest is exon end
	i <- which(isOnChrom & exons$extract(,"end") + 1L == position)
	if(length(i) > 0L) {
		anno <- rbind(
			anno,
			unique(with(exons$extract(i), data.frame(exon=groupPosition, end="right")))
		)
	}
	
	# Merge
	if(is.null(anno)) { exon <- NA
	} else            { exon <- with(anno[ order(anno$exon) ,],   paste(sprintf(ifelse(end == "right", "%s]", "[%s"), exon), collapse=","))
	}
	
	return(exon)
}

# Process one group of events (1 junction of interest + alternatives) of interest
processExtended <- function(x, exons) {
	# Junction annotation
	for(i in 1:nrow(x)) {
		# Exon boundary at splicing sites
		exon.left  <- annotateSplicingSite(chrom=x[i,"chrom"], position=x[i,"left"],  exons)
		exon.right <- annotateSplicingSite(chrom=x[i,"chrom"], position=x[i,"right"], exons)
		
		# site / partner rather than left / right
		if(x[i,"site.is.ID"] == "left") {
			x[i,"exon.site"] <- exon.left
			x[i,"exon.partner"] <- exon.right
		} else {
			x[i,"exon.site"] <- exon.right
			x[i,"exon.partner"] <- exon.left
		}
	}
	
	# Junction labelling (from most to least represented, except A = junction of interest)
	n <- apply(x[, grep("^I\\.", colnames(x)) ], 1, sum)
	x <- x[ order(n, decreasing=TRUE) ,]
	x[ x$ID != x$target , "label" ] <- LETTERS[ 1:sum(x$ID != x$target) + 1 ]
	x[ x$ID == x$target , "label" ] <- "A"
	
	# Sort
	x <- x[ order(x$site, x$ID, x$chrom, x$left, x$right) ,]
	
	# Samples with significant results on alternative A
	samples <- apply(x[ x$label == "A" , grep("^filter\\.", colnames(x)) ], 2, any)
	samples <- sub("^filter\\.", "", names(which(samples)))
	
	# Mark to plot gene in normalized coordinates
	symbol <- unique(x$ID.symbol)
	symbol <- symbol[ !is.na(symbol) & symbol != "" ]
	toPlot <- data.frame(
		sample = samples,
		symbol = symbol,
		stringsAsFactors = FALSE
	)
	
	return(list(x=x, toPlot=toPlot))
}

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
		dir.create(trackDir, showWarnings=FALSE)
		saveRDT(trk, file=trackFile)
	}
	
	return(trk)
}

# Plots junctions over a simplified representation of the transcript (normalized exon and intron sizes)
plot.normalized <- function(evt, sample, symbol, exons, outDir="out", bamDir="out/BAM", trackDir="out/depth", shape=1) {
	# BAM file
	bamFile <- sprintf("%s/%s.DNA.MD.sort.bam", bamDir, sample)
	if(!file.exists(bamFile)) stop("\"", bamFile, "\" doesn't exist")
	
	# Annotation of the transcript of interest
	gene <- exons$extract(exons$extract(,"symbol") == symbol)
	gene$ID <- paste(gene$start, gene$end, sep="-")
	
	# Exons (without redundancy)
	ano <- unique(gene[, c("start","end") ])
	
	# Genomic intervals
	breaks <- sort(unique(c(ano$start, ano$end)))
	ano <- data.frame(
		start = head(breaks, -1),
		end = tail(breaks, -1)
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
	x <- which(apply(!is.na(ano[,-c(1,2)]), 1, any))
	for(i in x) {
		trk <- depth(sample, bamFile, chrom=gene[1,"chrom"], start=ano[i,"start"], end=ano[i,"end"], trackDir=trackDir)
		ano[i,"depth"] <- sum(with(trk$extract(), (end-start+1L)*value)) / (ano[i,"end"] - ano[i,"start"] + 1L)
	}
	
	# Supporting reads
	evt$reads <- evt[, sprintf("I.%s", sample) ]
	
	# Filter status
	evt$filter <- evt[, sprintf("filter.%s", sample) ]
	
	# Restrict to supported junctions
	e <- evt[ evt$ID.symbol == symbol & evt$reads > 0L ,]
	
	# Normalized coordinates
	for(i in 1:nrow(e)) {
		for(site in c("left", "right")) {
			# Overlapping feature
			ovl <- which(ano$start < e[i,site] & ano$end >= e[i,site])
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
	
	# Image file
	width <- 100 + nrow(ano) * 30
	height <- 430 + length(transcripts) * 40
	file <- sprintf("%s/%s - %s.png", outDir, symbol, sample)
	png(file=file, width=width, height=height, res=100)

	# Layout
	layout(matrix(1:3, ncol=1), heights=c(lcm(4), 1, lcm(4)))
	par(oma=c(1,1,3,1), cex=1)
	xlim <- c(-0.5, nrow(ano)+0.5)
	
	# Annotated junctions
	if(any(e$class == "annotated")) { ymax <- log(max(e[ e$class == "annotated" , "reads" ]), 10)
	} else                          { ymax <- 1
	}
	par(mar=c(0,5,0,0))
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
	
	# Legend
	legend(
		x="top", horiz=TRUE, inset=-0.3, xpd=NA, bty="n",
		lwd = c(2, 2, 2, 2, 1),
		lty = c("solid", "solid", "solid", "solid", "dotted"),
		col = c("royalblue", "forestgreen", "orange", "red", "black"),
		legend = c("annotated", "plausible", "anchored", "unknown", "filtered")
	)
	
	# Transcripts
	par(mar=c(1,5,1,0))
	plot(x=NA, y=NA, xlim=xlim, ylim=c(0, length(transcripts)), xlab="", ylab="", xaxs="i", xaxt="n", yaxt="n", yaxs="i", bty="n")
	for(i in 1:length(transcripts)) {
		# Included exons
		x <- ano[[ transcripts[i] ]]
		
		# Plot
		level <- ano[ which(!is.na(x)) , "depth" ] / max(ano$depth, na.rm=TRUE)
		segments(x0=head(which(!is.na(x)), 1)-0.5, x1=tail(which(!is.na(x)), 1)-0.5, y0=i-0.5, y1=i-0.5)
		rect(xleft=which(!is.na(x))-1, xright=which(!is.na(x)), ybottom=i-0.9, ytop=i-0.1, col=grey(1 - level), border="#000000")
		text(x=which(!is.na(x))-0.5, y=i-0.5, adj=c(0.5, 0.5), labels=x[!is.na(x)], col=ifelse(level < 0.5, "black", "white"))
		mtext(side=2, at=i-0.5, text=sub(" .+$", "", transcripts[i]), las=2, line=0)
	}
	
	# Other junctions
	if(any(e$class != "annotated")) { yaxt="s"; ymax <- log(max(e[ e$class != "annotated" , "reads" ]), 10)
	} else                          { yaxt="n"; ymax <- 1
	}
	par(mar=c(0,5,0,0))
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
	
	dev.off()

}

# Export detailed table in CSV
exportDetails <- function(tab, file="out/Details.csv") {
	# Column selection
	exp <- tab[, c("target", "site", "label", "chrom", "left", "right", "ID.symbol", "class", "exon.site", "exon.partner") ]
	colnames(exp) <- c("Jonction d'intérêt", "Site d'épissage", "Alternative", "Chrom", "Gauche", "Droite", "Gene", "Classe", "Exon (site)", "Exon (partenaire)")
	
	# Rename I and round PSI
	for(sample in sub("^PSI\\.", "", grep("^PSI", colnames(tab), value=TRUE))) {
		exp[[ sprintf("%s.I", sample) ]] <- tab[[ sprintf("I.%s", sample) ]]
		exp[[ sprintf("%s.PSI", sample) ]] <- round(tab[[ sprintf("PSI.%s", sample) ]], 3)
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
	
	# Print data
	startRow <- 1L
	for(target in unique(details$"Jonction d'intérêt")) {
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
			tmp <- tab[ rows , c("target", "site", "chrom", "left", "right", "ID.symbol", "class", "exon.site", "exon.partner", sprintf("PSI.%s", sample), sprintf("I.%s", sample)) ]
			tmp$is.left <- tmp$site == sprintf("%s:%i", tmp$chrom, tmp$left)
			tmp$site <- NULL
			
			# Merge row pairs
			mrg <- tmp[ tmp$is.left ,]
			mrg$is.left <- NULL
			mrg$sample <- sample
			mrg$reads <- mrg[[ sprintf("I.%s", sample) ]]
			mrg$PSI.left <- tmp[  tmp$is.left , sprintf("PSI.%s", sample) ]
			mrg$PSI.right <- tmp[ !tmp$is.left , sprintf("PSI.%s", sample) ]
			mrg[[ sprintf("PSI.%s", sample) ]] <- NULL
			mrg[[ sprintf("I.%s", sample) ]] <- NULL
			colnames(mrg)[ colnames(mrg) == "exon.site" ] <- "exon.left"
			colnames(mrg)[ colnames(mrg) == "exon.partner" ] <- "exon.right"
			
			# Aggregate
			cand[[ sample ]] <- mrg
		}
	}
	
	# Merge samples
	cand <- do.call(rbind, cand)
	
	# Prioritize
	cand <- cand[ order(cand$reads, decreasing=TRUE) ,]
	
	# Rename columns
	cand <- cand[, c("target", "chrom", "left", "right", "ID.symbol", "class", "exon.left", "exon.right", "sample", "reads", "PSI.left", "PSI.right") ]
	colnames(cand) <- c("Jonction d'intérêt", "Chrom", "Gauche", "Droite", "Gene", "Classe", "Exon (gauche)", "Exon (droite)", "Echantillon", "Reads", "PSI (gauche)", "PSI (droite)")
	
	# Export
	if(!is.na(file)) write.csv2(cand, file=file, row.names=FALSE, na="")
	
	invisible(cand)
}



message("Loading dependencies...")

library(Rgb)
library(openxlsx)

message("Parsing annotation...")

exons <- readRDT(exonFile)

message("Parsing events...")

events <- readRDS(eventFile)

message("Filtering events...")

events <- filterEvents(events, min.PSI=min.PSI, min.I=min.I)

message("Extending events...")

extended <- extendEvents(events)

message("Filtering event groups...")

extended.keep <- filterExtended(extended, symbols.filter=symbolList, A.types.filter=classList, focus=focusList)
message("- ", sum(extended.keep), " junctions of interest")

message("Processing event groups...")

out <- list()
toPlot <- list()
for(junction in names(which(extended.keep))) {
	
	message("- ", junction)
	
	# Subset site
	x <- extended[[ junction ]]
	x$target <- junction
	x$site <- factor(x$site)
	
	# Process
	tmp <- processExtended(x, exons=exons)
	out[[ length(out) + 1L ]] <- tmp$x
	toPlot[[ length(toPlot) + 1L ]] <- tmp$toPlot
}

message("Merging output...")

out <- do.call(rbind, out)

message("Preparing output directory...")

outDir <- sprintf("I-%i_PSI-%g_%s_%s_%s", min.I, min.PSI, symbols, classes, gsub(":", "-", focus))
dir.create(outDir)

if(isTRUE(plot) && length(toPlot) > 0L) {
	
	message("Merging genes to plot...")
	
	toPlot <- do.call(rbind, toPlot)
	toPlot <- unique(toPlot)
	
	message("Plotting ", nrow(toPlot), " genes & samples ...")
	
	dir.create(sprintf("%s/plots", outDir))
	for(i in 1:nrow(toPlot)) {
		message("- ", toPlot[i,"symbol"], " - ", toPlot[i,"sample"])
		plot.normalized(evt=events, sample=toPlot[i,"sample"], symbol=toPlot[i,"symbol"], exons=exons, outDir=sprintf("%s/plots", outDir), bamDir=".", trackDir="depth")
	}
} else message("Nothing to plot")

message("Exporting candidates...")

exportCandidates(out, file=sprintf("%s/Candidates.csv", outDir))

message("Exporting details...")

details <- exportDetails(out, file=sprintf("%s/Details.csv", outDir))
formatDetails(details, out, file=sprintf("%s/Details.xlsx", outDir))

message("done")
