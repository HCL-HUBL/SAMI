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
  "            - UMI.median:",
  "                namespace: 'UMI'",
  "                description: 'Median UMI family size (reads sharing similar UMI and genomic coordinates)'",
  "                format: '{:,.1f}'",
  "            - UMI.unique:",
  "                namespace: 'UMI'",
  "                description: 'Proportion of reads which were not duplicated (unique combinations of UMI and genomic coordinates)'",
  "                suffix: '%'",
  "                format: '{:,.2f}'",
  "        data:"
)
cat(lines, sep="\n", file="./umi_table.yaml")

### Data
for(ifile in allfile)
{
  samp <- gsub(x=ifile, pattern="_family_size_histogram.txt", replacement="")
  umi_hist <- read.delim(file=ifile)
  tostat <- rep(x=umi_hist$family_size, times=umi_hist$count)
  ### Compute the "true" proportion of UMI with one read
  trueOne <- 100 * (umi_hist$count*umi_hist$family_size / sum(umi_hist$count*umi_hist$family_size))[1]
  lines <- sprintf("            '%s':
                UMI.median: %f
                UMI.unique: %f", samp, median(tostat), trueOne)
            
  cat(lines, sep="\n", file="./umi_table.yaml", append=TRUE)
}
