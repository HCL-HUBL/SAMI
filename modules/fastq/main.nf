process fastq {
	tag "$sample"

	cpus 1
	time { 5.minute * task.attempt }
	memory { 500.MB * task.attempt }

	input:
	tuple path(R1), path(R2), val(sample), val(pair), val(type)
	val(CN)
	val(PL)
	val(PM)

	output:
	tuple path(R1), path(R2), val(sample), val("${type.first()}"), stdout, emit: FASTQ

	"""
	#!/usr/bin/env Rscript --vanilla

	# Get FASTQ sets from Nextflow
	R1 <- strsplit("${R1.join("|")}", split="|", fixed=TRUE)[[1]]
	R2 <- strsplit("${R2.join("|")}", split="|", fixed=TRUE)[[1]]

	# Check type consistency
	type <- strsplit("${type.join("|")}", split="|", fixed=TRUE)[[1]]
	type <- unique(type)
	if(length(type) != 1L) stop("A sample can't mix single and paired-end FASTQ")
	pairedEnd <- type == "paired"

	# Parse regex list
	def <- scan("${moduleDir}/etc/FASTQ_headers.txt", what="", sep="\n", quiet=TRUE)
	regex <- def[ c(TRUE,FALSE) ]
	names <- strsplit(def[ c(FALSE,TRUE) ], split=" ")

	# For each R1/R2 pair
	RG <- NULL
	for(i in 1:length(R1)) {

	# Get first read headers (whether the file is compressed or not)
	H1 <- scan(R1[i], what="", sep="\n", n=1L, quiet=TRUE)
	H2 <- scan(R2[i], what="", sep="\n", n=1L, quiet=TRUE)
	if(length(H1) == 0L)              stop("No header in R1")
	if(length(H2) == 0L && pairedEnd) stop("No header in R2")

	# For each defined regex
	metadata_R1 <- character(0)
	metadata_R2 <- character(0)
	for(j in 1:length(regex)) {
		# Matching regex (R1)
		if(length(metadata_R1) == 0L && grepl(regex[j], H1)) {
		# Extract elements
		metadata_R1 <- regmatches(H1, regexec(regex[j], H1))[[1]][-1]
		names(metadata_R1) <- names[[j]]
		}

		# Matching regex (R2)
		if(pairedEnd && length(metadata_R2) == 0L && grepl(regex[j], H2)) {
		# Extract elements
		metadata_R2 <- regmatches(H2, regexec(regex[j], H2))[[1]][-1]
		names(metadata_R2) <- names[[j]]
		}
	}

	# No match
	if(length(metadata_R1) == 0L)              stop("Unable to parse the header of R1")
	if(pairedEnd && length(metadata_R2) == 0L) stop("Unable to parse the header of R2")

	# Run identity
	ID1 <- metadata_R1[ c("instrument", "run", "flowcell", "lane", "index") ]
	if(pairedEnd) {
	ID2 <- metadata_R2[ c("instrument", "run", "flowcell", "lane", "index") ]
	if(!identical(ID1, ID2)) stop("Wrong R1/R2 FASTQ file matching")
	}

	# Mate identity (if any is provided, both are required)
	R1_has_read_info <- !is.na(metadata_R1["read"]) && metadata_R1["read"] != ""
	R2_has_read_info <- !is.na(metadata_R2["read"]) && metadata_R2["read"] != ""
	if(R1_has_read_info || R2_has_read_info) {
	if(is.na(metadata_R1["read"]) || metadata_R1["read"] != "1")                   stop("R1 FASTQ file does not contain R1 reads")
	if(pairedEnd && (is.na(metadata_R2["read"]) || !metadata_R2["read"] %in% 2:3)) stop("R2 FASTQ file does not contain R2 or R3 reads")
	}

	# BC is optional
	if(is.na(metadata_R1["index"]) || metadata_R1["index"] == "") { BC <- ""
	} else                                                        { BC <- sprintf("\tBC:%s", metadata_R1["index"])
	}

	# RG definition for BAM header (current pair)
	x <- sprintf("ID:%s_%i%s", "${sample}", i, BC)
	if("${CN}" != "") x <- c(x, "CN:${CN}")
	if("${PL}" != "") x <- c(x, "PL:${PL}")
	if("${PM}" != "") x <- c(x, "PM:${PM}")
	x <- c(x, sprintf("PU:%s", paste(ID1, collapse=":")))
	x <- c(x, "SM:${sample}")

	# Merge with all read pairs
	RG <- c(RG, paste(x, collapse="\t"))
	}

	# Print final RG to stdout
	cat(paste(RG, collapse=" , "))
	"""
}
