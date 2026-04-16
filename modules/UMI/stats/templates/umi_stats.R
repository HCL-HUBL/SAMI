#!/usr/bin/env Rscript --vanilla

# Source files
files <- dir("./histograms", pattern=".*_family_size_histogram.txt", full.names=TRUE)

# General stats - header
gstats <- c(
	"custom_data:",
	"    umi_stats:",
	"        plot_type: 'generalstats'",
	"        headers:",
	"            - UMI.median:",
	"                namespace: 'UMI'",
	"                description: 'Median UMI family size (reads sharing similar UMI and genomic coordinates)'",
	"                format: '{:,.1f}'",
	"            - UMI.unique:",
	"                namespace: 'UMI'",
	"                description: 'Proportion of reads which were never duplicated (unique combinations of UMI and genomic coordinates)'",
	"                suffix: '%'",
	"                format: '{:,.2f}'",
	"            - UMI.dup:",
	"                namespace: 'UMI.dup'",
	"                description: 'Proportion of reads considered as duplicates of other reads'",
	"                format: '{:,.2f}'",
	"                min: 0",
	"                max: 100",
	"                suffix: '%'",
	"        data:"
)

# Plot - header
plot <- c(
	"id: 'UMI_duplication'",
	"section_name: 'UMI duplication'",
	"description: 'distribution of the number of UMIs depending on their number of copies.'",
	"plot_type: 'linegraph'",
	"headers:",
	"    id: 'UMI_duplication_linegraph'",
	"    title: 'UMI duplication'",
	"    xlab: 'Number of UMI copies'",
	"    ylab: 'Fraction of the total number of read'",
	"data:"
)

for(file in files) {
	# Parse histogram
	sample <- sub("_family_size_histogram.txt$", "", basename(file))
	histo <- read.delim(file)
	total <- sum(histo$count * histo$family_size)
	
	# UMI.median
	x <- rep(x=histo$family_size, times=histo$count)
	UMI.median <- median(x)
	
	# UMI.unique
	UMI.unique <- 100 * histo[ histo$family_size == 1L , "count" ] / total
	
	# UMI.dup
	UMI.dup <- 100 - 100 * sum(histo$count) / total
	
	# Fraction
	fraction <- histo$count * histo$family_size / total
	
	# General stats - content
	gstats <- c(
		gstats,
		sprintf("            '%s':", sample),
		sprintf("                UMI.median: %f", UMI.median),
		sprintf("                UMI.unique: %f", UMI.unique),
		sprintf("                UMI.dup: %f", UMI.dup)
	)
	
	# Plot - content
	values <- paste(
		sprintf(
			"%i: %g",
			histo$family_size,
			fraction
		),
		collapse = ", "
	)
	plot <- c(
		plot,
		sprintf(
			"    %s: { %s }",
			sample,
			values
		)
	)
}

# General stats - write
cat(gstats, sep="\n", file="./umi-table.yaml")

# Plot - write
cat(plot, sep="\n", file="./umi-plot_mqc.yaml")
