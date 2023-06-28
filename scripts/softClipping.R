#!/usr/bin/env Rscript --vanilla

# CLI arguments
args <- commandArgs(TRUE)
if(length(args) != 2L) stop("USAGE: ./softClipping.R SAMPLE_ID BAM_FILE")
sample <- args[1]
bamFile <- args[2]

# Check arguments
if(!file.exists(bamFile)) stop("BAM_FILE must exist")



# Collect CIGARs (R1 = forward, R2 = reverse)
r1 <- scan(pipe(sprintf("samtools view -F 0x10 \"%s\" | cut -f6", bamFile)), what="", sep="\n", quiet=TRUE)
r2 <- scan(pipe(sprintf("samtools view -f 0x10 \"%s\" | cut -f6", bamFile)), what="", sep="\n", quiet=TRUE)

# Count soft-clipped bases on each side
r1.left <- as.integer(sub("^([0-9]+)S.*$", "\\1", r1))
r2.left <- as.integer(sub("^([0-9]+)S.*$", "\\1", r2))
r1.right <- as.integer(sub("^.*?([0-9]+)S$", "\\1", r1))
r2.right <- as.integer(sub("^.*?([0-9]+)S$", "\\1", r2))

# NA is no-soft-clipping
r1.left[ is.na(r1.left) ] <- 0L
r2.left[ is.na(r2.left) ] <- 0L
r1.right[ is.na(r1.right) ] <- 0L
r2.right[ is.na(r2.right) ] <- 0L

# Summary
smax <- max(r1.left, r2.left, r1.right, r2.right)
out <- rbind(
	R1.5p = tabulate(r1.left+1L, nbins=smax),
	R1.3p = tabulate(r1.right+1L, nbins=smax),
	R2.5p = tabulate(r2.right+1L, nbins=smax),
	R2.3p = tabulate(r2.left+1L, nbins=smax)
)
colnames(out) <- 1:ncol(out) - 1L

# Normalize
out <- out / sum(out[1,])

x <- as.integer(colnames(out))
for(read in c("R1", "R2")) {
	for(side in c("5p", "3p")) {
		# YAML graphs
		y <- round(100*out[ sprintf("%s.%s", read, side) ,], 3)
		lines <- c(
			"parent_id: 'soft_clipping'",
			"parent_name: 'Soft-clipping'",
			"parent_description: 'Soft-clipping observed during alignment. Constant values suggest the presence of untrimmed adapters or UMIs.'",
			sprintf("id: 'Soft_clipping_%s_%s'", read, side),
			sprintf("section_name: '%s %s'", read, side),
			sprintf("description: 'Soft-clipping occuring on the %s side of read %s during alignment.'", side, read),
			"plot_type: 'linegraph'",
			"pconfig:",
			sprintf("    id: 'Soft_clipping_%s_%s_linegraph'", read, side),
			sprintf("    title: 'Soft clipping (%s, %s)'", read, side),
			"    xlab: 'Soft clipping length (bp)'",
			"    ylab: 'Proportion of all aligned reads (%)'",
			"data:",
			sprintf("    '%s': { %s }", sample, paste(sprintf("%i: %g", x, y), collapse=", "))
		)
		cat(lines, sep="\n", file=sprintf("%s_%s-%s_mqc.yaml", sample, read, side))
	}
}

