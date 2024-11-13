#!/usr/bin/env Rscript --vanilla

### Read the list of files
allfile <- list.files(path="./", pattern=".*_family_size_histogram.txt")

### YAML table
### Header
lines <- c(
  "custom_data:",
  "    umi_stats:",
  "        plot_type: 'generalstats'",
  "        headers:",
  "            - UMI.meanBeforeConsensus:",
  "                namespace: 'UMI'",
  "                description: 'Mean, before consensus, of the number of UMIs per sequence'",
  "                format: '{:,.2f}'",
  "            - UMI.medianBeforeConsensus:",
  "                namespace: 'UMI'",
  "                description: 'Median, before consensus, of the number of UMIs per sequence'",
  "                format: '{:,.1f}'",
  "            - UMI.meanAfterConsensus:",
  "                namespace: 'UMI'",
  "                description: 'Mean, after consensus, of the number of UMIs per sequence'",
  "                format: '{:,.2f}'",
  "            - UMI.medianAfterConsensus:",
  "                namespace: 'UMI'",
  "                description: 'Median, after consensus, of the number of UMIs per sequence'",
  "                format: '{:,.1f}'",
  "            - UMI.max:",
  "                namespace: 'UMI'",
  "                description: 'Maximum number of UMIs per sequence'",
  "                format: '{:,.0f}'",
  "            - UMI.unique:",
  "                namespace: 'UMI'",
  "                description: 'Percentage of sequences with one read per UMI'",
  "                suffix: '%'",
  "                format: '{:,.2f}'",
  "        data:"
)
cat(lines, sep="\n", file="./umi_table_mqc.yaml")

### Data
for(ifile in allfile)
{
  samp <- gsub(x=ifile, pattern="_family_size_histogram.txt", replacement="")
  umi_hist <- read.delim(file=ifile)
  tostat <- rep(x=umi_hist$family_size, times=umi_hist$count)
  tostat2 <- rep(x=umi_hist$family_size, times=(umi_hist$count*umi_hist$family_size))
  ### Compute the "true" proportion of UMI with one read
  trueOne <- 100 * (umi_hist$count*umi_hist$family_size / sum(umi_hist$count*umi_hist$family_size))[1]
  lines <- sprintf("            '%s':
                UMI.meanBeforeConsensus: %f
                UMI.medianBeforeConsensus: %f
                UMI.meanAfterConsensus: %f
                UMI.medianAfterConsensus: %f
                UMI.max: %d
                UMI.unique: %f", samp, mean(tostat2), median(tostat2), mean(tostat), median(tostat), max(tostat), trueOne)
            
  cat(lines, sep="\n", file="./umi_table_mqc.yaml", append=TRUE)
}
