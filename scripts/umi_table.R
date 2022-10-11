### Read the list of files
allfile <- list.files(path="./", pattern=".*_family_size_histogram.txt")

### YAML table
### Header
lines <- c(
  "custom_data:",
  "    umi_stats:",
  "        plot_type: 'generalstats'",
  "        pconfig:",
  "            - UMI.mean:",
  "                namespace: 'UMI'",
  "                description: 'Mean number of UMIs per sequence'",
  "                format: '{:,.0f}'",
  "            - UMI.median:",
  "                namespace: 'UMI'",
  "                description: 'Median number of UMIs per sequence'",
  "                format: '{:,.0f}'",
  "            - UMI.max:",
  "                namespace: 'UMI'",
  "                description: 'Maximum number of UMIs per sequence'",
  "                format: '{:,.0f}'",
  "            - UMI.unique:",
  "                namespace: 'UMI'",
  "                description: 'Proportion of sequences with one read per UMI'",
  "                format: '{:,.0f}'",
  "        data:"
)
cat(lines, sep="\n", file="./umi_table_mqc.yaml")

### Data
for(ifile in allfile)
{
  samp <- gsub(x=ifile, pattern="_family_size_histogram.txt", replacement="")
  umi_hist <- read.delim(file=ifile)
  tostat <- rep(x=umi_hist$family_size, times=umi_hist$count)
  lines <- c(sprintf("            '%s':
                UMI.mean: %f
                UMI.median: %f
                UMI.max: %f,
                UMI.unique: %f", samp, mean(tostat), median(tostat), max(tostat), umi_hist$fraction[1])
)
  cat(lines, sep="\n", file="./umi_table_mqc.yaml", append=TRUE)
}
