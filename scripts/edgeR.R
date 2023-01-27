#!/usr/bin/env Rscript --vanilla

### Collect individual featureCount outputs and compute QC using edgeR
### Author : <sylvain.mareschal@lysarc.org>

# Collect CLI arguments
args <- commandArgs(TRUE)
if(length(args) < 3L) stop("USAGE : ./edgeR.R ANNOTATION.csv OUT_DIR COUNT_FILE_1 [ COUNT_FILE_2 [ ... ] ]")
annotation <- args[1]
outDir <- args[2]
countFiles <- args[ 3:length(args) ]

# Check CLI arguments
if(!file.exists(annotation))      stop("ANNOTATION.csv must exist")
if(!file.exists(outDir))          stop("OUT_DIR must exist")
if(any(!file.exists(countFiles))) stop("COUNT_FILE_N must all exist")

# Dependencies
library(edgeR)



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
write.csv2(counts, file=sprintf("%s/counts.csv", outDir), row.names=TRUE, col.names=NA)



# Gene size (for RPKM)
sizes <- read.csv(annotation, stringsAsFactors=FALSE, row.names="GeneID")
sizes <- sizes[ rownames(counts) , "Length" ]
names(sizes) <- rownames(counts)

# Import into EdgeR
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge, method="TMM")

# Compute RPKs for all genes
RPK <- counts / (sizes / 1e3)

# Export CPM matrix
write.csv2(RPK, file=sprintf("%s/RPK.csv", outDir), row.names=TRUE, col.names=NA)

# Compute (normalized) CPMs for genes of interest
cpm <- cpm(dge)

# Export CPM matrix
write.csv2(cpm, file=sprintf("%s/CPM.csv", outDir), row.names=TRUE, col.names=NA)

# Classify RPKs
class <- matrix(cut(RPK, breaks=c(0L, 1L, 5L, 10L, 20L, 100L)), nrow=nrow(RPK), ncol=ncol(RPK), dimnames=dimnames(RPK))
class[ counts == 0L ] <- "No read"
class[ RPK > 100L ] <- ">100"

# Summarize RPKs
mtx <- NULL
for(i in 1:ncol(class)) mtx <- rbind(mtx, table(factor(class[,i], levels=c("No read", "(0,1]", "(1,5]", "(5,10]", "(10,20]", "(20,100]", ">100"))))
rownames(mtx) <- colnames(class)

# Reliable genes
## VW: modification "drop=FALSE" to take into account the 1 sample case
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
	"        pconfig:",
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
	"pconfig:",
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

