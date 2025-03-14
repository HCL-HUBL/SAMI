#!/usr/bin/env Rscript --vanilla

### Collect Nextflow arguments
sample    <- "!{sample}"
inputFile <- "!{umiHist}"

if(!file.exists(inputFile)) stop("INPUT_FILE must exist")

### Read the file
umi_hist <- read.delim(file=inputFile)

###
ori.fraction <- umi_hist$count * umi_hist$family_size / sum(umi_hist$count * umi_hist$family_size)

### YAML graph
lines <- c(
  "id: 'UMI_duplication'",
  "section_name: 'UMI duplication'",
  "description: 'distribution of the number of UMIs depending on their number of copies.'",
  "plot_type: 'linegraph'",
  "headers:",
  "    id: 'UMI_duplication_linegraph'",
  "    title: 'UMI duplication'",
  "    xlab: 'Number of UMI copies'",
  "    ylab: 'Fraction of the total number of read'",
  "data:",
  sprintf("    %s: { %s }", sample, paste(sprintf("%i: %g", umi_hist$family_size, ori.fraction), collapse=", "))
)
cat(lines, sep="\n", file=sprintf("./%s_mqc.yaml", sample))
