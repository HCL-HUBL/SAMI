#!/usr/bin/env Rscript --vanilla

# Collect Nextflow arguments
range <- as.integer("!{range}")
exonFile <- "!{exons}"
fastaFile <- "!{fasta}"
indexFile <- "!{fai}"



library(Rgb)

mrg <- NULL
files <- dir("gapFiles", pattern="\\.out\\.tab$", full.names=TRUE)
for(file in files) {
	# Parse file
	tab <- read.table(
		file, sep="\t", quote=NULL, comment.char="",
		col.names = c("chrom", "start", "end", "strand", "motif", "annotated", "reads.uni", "reads.multi", "overhang"),
		colClasses = c("character", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer"),
	)
	
	# Junction ID
	rownames(tab) <- sprintf("%s:%i-%i", tab$chrom, tab$start, tab$end)
	
	# Merge
	if(is.null(mrg)) {
		# Start with first table
		mrg <- tab
	} else {
		# Add read counts for known junctions
		common <- intersect(rownames(tab), rownames(mrg))
		mrg[ common , "reads.uni" ] <- mrg[ common , "reads.uni" ] + tab[ common , "reads.uni" ]
		mrg[ common , "reads.multi" ] <- mrg[ common , "reads.multi" ] + tab[ common , "reads.multi" ]
		mrg[ common , "overhang" ] <- max(mrg[ common , "overhang" ], tab[ common , "overhang" ])
		
		# Add novel junctions
		new <- setdiff(rownames(tab), rownames(mrg))
		mrg <- rbind(mrg, tab[new,])
	}
}

# Resort genomically
mrg <- mrg[ order(mrg$chrom, mrg$start, mrg$end) ,]

# Exon annotation
exons <- readRDT(exonFile)

# Reference genome
fasta <- track.fasta.multi(fastaFile, indexFile)

# Exon boundaries
exons.left <- unique(exons$extract(,c("chrom", "start")))
exons.left <- split(exons.left[[2]], exons.left[[1]])
exons.right <- unique(exons$extract(,c("chrom", "end")))
exons.right <- split(exons.right[[2]], exons.right[[1]])

# Unannotated junctions to fix
indexes <- which(mrg$annotated == 0L)
for(index in indexes) {
	# Chromosome
	chrom <- sub("^chr", "", mrg[index,"chrom"])
	
	# Is there an exon end near the gap start ?
	shift.left <- exons.right[[ chrom ]] - mrg[index,"start"] + 1L
	if(any(shift.left == 0L)) next;
	exon.left <- which(abs(shift.left) <= range)
	
	# Is there an exon start near the gap end ?
	shift.right <- exons.left[[ chrom ]] - mrg[index,"end"] - 1L
	if(any(shift.right == 0L)) next;
	exon.right <- which(abs(shift.right) <= range)
	
	# Candidate shifts to assess
	shifts <- unique(c(
		shift.left[exon.left],
		shift.right[exon.right]
	))
	
	for(shift in shifts) {
		if(shift > 0L) {
			# First bases of the gap
			left <- fasta$slice(chrom=mrg[index,"chrom"], start=mrg[index,"start"], end=mrg[index,"start"] + shift - 1L, multiple=FALSE)
			
			# Bases after the gap
			right <- fasta$slice(chrom=mrg[index,"chrom"], start=mrg[index,"end"] + 1L, end=mrg[index,"end"] + shift, multiple=FALSE)
		} else {
			# Bases before the gap
			left <- fasta$slice(chrom=mrg[index,"chrom"], start=mrg[index,"start"] + shift, end=mrg[index,"start"] - 1L, multiple=FALSE)
			
			# Last bases of the gap
			right <- fasta$slice(chrom=mrg[index,"chrom"], start=mrg[index,"end"] + shift + 1L, end=mrg[index,"end"], multiple=FALSE)
		}
		
		if(left == right) {
			# Apply shift
			mrg[index,"start"] <- mrg[index,"start"] + shift
			mrg[index,"end"] <- mrg[index,"end"] + shift
			
			message(rownames(mrg)[index], " : shift ", shift, " OK")
		} else {
			message(rownames(mrg)[index], " : shift ", shift, " ERROR")
		}
	}
}

# Export
write.table(mrg, file="fixed.out.tab", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
