#!/usr/bin/env Rscript --vanilla

# Collect Nextflow arguments
annotation <- "!{annotation}"  
countFiles <- strsplit("!{countFiles.join('|')}", split="|", fixed=TRUE)[[1]]

# Dependencies
library(edgeR)



# Output directory
outDir <- "."

# Aggregate a count matrix
counts <- NULL
for(countFile in countFiles) {
	# Parse counts
	tmp <- readRDS(countFile)
	sample <- colnames(tmp)
	message(sample)
	
	# Aggregate
	if(is.null(counts)) { counts <- tmp
	} else              { counts <- cbind(counts, tmp)
	}
}

# Export count matrix
write.table(counts, file=sprintf("%s/counts.tsv", outDir), row.names=TRUE, col.names=NA, sep="\t")



# Gene size (for RPKM)
sizes <- read.table(annotation, stringsAsFactors=FALSE, header=TRUE, row.names="GeneID", sep="\t")
sizes <- sizes[ rownames(counts) , "Length" ]
names(sizes) <- rownames(counts)

# Import into EdgeR
dge <- DGEList(counts=counts)

# Compute normalization factors (for samples with counts)
empty <- apply(counts, 2L, sum) == 0L
factors <- rep(as.numeric(NA), ncol(dge))
factors[!empty] <- calcNormFactors(dge$counts[,!empty], method="TMM")
dge$samples$norm.factors <- factors

# Compute RPKs for all genes
RPK <- counts / (sizes / 1e3)

# Export CPM matrix
write.table(RPK, file=sprintf("%s/RPK.tsv", outDir), row.names=TRUE, col.names=NA, sep="\t")

# Compute (normalized) CPMs for genes of interest (for samples with counts)
cpm <- cpm(dge[,!empty])
if(any(empty)) {
       extra <- matrix(0L, nrow=nrow(cpm), ncol=sum(empty), dimnames=list(rownames(cpm), colnames(dge)[empty]))
       cpm <- cbind(cpm, extra)
       cpm <- cpm[, colnames(dge) ]
}

# Export CPM matrix
write.table(cpm, file=sprintf("%s/CPM.tsv", outDir), row.names=TRUE, col.names=NA, sep="\t")

# Classify RPKs
class <- matrix(cut(RPK, breaks=c(0L, 1L, 5L, 10L, 20L, 100L)), nrow=nrow(RPK), ncol=ncol(RPK), dimnames=dimnames(RPK))
class[ counts == 0L ] <- "No read"
class[ RPK > 100L ] <- ">100"

# Summarize RPKs
mtx <- NULL
for(i in 1:ncol(class)) mtx <- rbind(mtx, table(factor(class[,i], levels=c("No read", "(0,1]", "(1,5]", "(5,10]", "(10,20]", "(20,100]", ">100"))))
rownames(mtx) <- colnames(class)

# Reliable genes
RPK5 <- apply(mtx[, c("(5,10]", "(10,20]", "(20,100]", ">100"), drop=FALSE ], 1, sum)

# Add edgeR library info
tab <- as.data.frame(mtx)
tab$lib.size <- dge$samples[ rownames(tab) , "lib.size" ]
norm.factors <- dge$samples[ rownames(tab) , "norm.factors" ]
tab$TMM.sum <- tab$lib.size * norm.factors
tab$Genes_RPK5 <- RPK5



### EXPORT GENERAL STATS ###

# Header     
lines <- c(
	"custom_data:",
	"    edgeR_stats:",
	"        plot_type: 'generalstats'",
	"        headers:",
	"            - lib.size:",
	"                namespace: 'edgeR'",
	"                description: 'Total counts (reads non-ambigously assigned to a transcript)'",
	"                format: '{:,.0f}'",
	"            - TMM.sum:",
	"                namespace: 'edgeR'",
	"                description: 'TMM-normalized library size'",
	"                format: '{:,.0f}'",
	"            - Genes_RPK5:",
	"                title: 'Reliable genes'",
	"                namespace: 'Gene_coverage'",
	"                description: 'Amount of genes with reliable gene expression estimates (at least 5 Reads Per Kilobase of transcript)'",
	"                format: '{:,.0f}'",
	"        data:"
)

# Data
tmp <- tab[, c("lib.size", "TMM.sum", "Genes_RPK5") ]
for(i in 1:nrow(tab)) {
	lines <- c(
		lines,
		sprintf("            '%s':", rownames(tmp)[i]),
		sprintf("                %s: %g", colnames(tmp), tmp[i,])
	)
}

# Export as YAML
cat(lines, sep="\n", file=sprintf("%s/edgeR.yaml", outDir))



### EXPORT SEQUENCING DEPTH (BARPLOT) ###

# Sort
tmp <- tab[ order(tab$"No read", decreasing=TRUE) ,]

# Write section file
lines <- c(
	"id: 'Gene_coverage'",
	"section_name: 'Gene coverage'",
	"description: 'assessment, in Reads Per Kilobase of transcript (RPK). A high proportion of lowly covered genes can be attributed to insufficient sequencing depth and could be rescued by resequencing the same library. On the other hand a suspiciously high proportion of medium-covered genes can be due to a poor quality library (random reads all along the genome). Unfortunately the optimum is unknown and sample-specific, keep in mind though that not all genes are expected to be expressed at a significant level in all conditions.'",
	"plot_type: 'bargraph'",
	"headers:",
	"    id: 'Sequencing_depth_bargraph'",
	"    title: 'Gene coverage (RPK)'",
	"data:",
	sprintf(
		"    '%s': {'No overlapping read (0)': %i, ']0;1]': %i, ']1;5]': %i, ']5;10]': %i, ']10;20]': %i, ']20;100]': %i, '>100': %i}",
		rownames(tmp),
		tmp$"No read",
		tmp$"(0,1]",
		tmp$"(1,5]",
		tmp$"(5,10]",
		tmp$"(10,20]",
		tmp$"(20,100]",
		tmp$">100"
	)
)

# Export as YAML
cat(lines, sep="\n", file=sprintf("%s/edgeR_mqc.yaml", outDir))

