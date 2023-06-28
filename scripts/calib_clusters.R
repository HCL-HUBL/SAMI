#!/usr/bin/env Rscript --vanilla

# CLI arguments
args <- commandArgs(TRUE)
if(length(args) != 2L) stop("USAGE: ./calib_clusters.R SAMPLE_ID FILE.fastq.gz")
sample <- args[1]
fastqFile <- args[2]

# Check arguments
if(!file.exists(fastqFile)) stop("FASTQ_FILE must exist")



# Amount of reads to process simultaneously
chunkSize <- 100000L

# Open gzipped FASTQ file
con <- gzfile(fastqFile, "r")
on.exit(close(con))

# Expected quality values
dic <- c("."=1L, "$"=2L, "A"=3L, "K"=4L)

# Reusable storage
class <- rep(as.character(NA), chunkSize)

# Output storage
levels <- c("All bases > 90%", "1 base < 90%", "2-5 bases < 90%", "Other", "Any base < 50%")
out <- rep(0L, length(levels))
names(out) <- levels

last <- FALSE
while(!last) {
	# Parse a chunk of reads
	fastq <- scan(con, what="", sep="\n", n=4*chunkSize, quiet=TRUE)
	
	# Focus on quality string
	qual <- fastq[ c(FALSE,FALSE,FALSE,TRUE) ]
	qual <- strsplit(qual, split="")
	
	# Last chunk
	if(length(qual) < chunkSize) {
		class <- class[1:length(qual)]
		last <- TRUE
	}
	
	# Count non-K positions
	nonK <- sapply(qual, function(x) { length(x) - sum(dic[x] == 4L) })
	
	# Classify clusters
	class[ nonK == 0L ] <- "All bases > 90%"
	class[ nonK == 1L ] <- "1 base < 90%"
	class[ nonK %in% 2:5 ] <- "2-5 bases < 90%"

	# Count . positions
	dot <- sapply(qual, function(x) { sum(dic[x] == 1L) })
	class[ dot > 0L ] <- "Any base < 50%"

	# Remainder
	class[ is.na(class) ] <- "Other"
	
	# Count classes
	tab <- table(class)
	out[ names(tab) ] <- out[ names(tab) ] + tab
}

lines <- c(
	"id: 'Calib_section'",
	"section_name: 'Calib'",
	"description: 'UMI cluster homogeneity. For each cluster of reads considered by Calib to originate from the same molecule, this QC counts the amount of positions with poor consensus. All clusters falling in \"All bases > 90%\" (over-clustering) or too many clusters with poor consensus (under-clustering) would advocate for better tuning of Calib parameters.'",
	"plot_type: 'bargraph'",
	"pconfig:",
	"    id: 'Calib_bargraph'",
	"    title: 'Calib clusters homogeneity'",
	"data:",
	sprintf("    %s: {%s}", sample, paste(sprintf("'%s': %i", names(out), out), collapse=", "))
)
cat(lines, sep="\n")

