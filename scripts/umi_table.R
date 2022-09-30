### Read the list of files
allfile <- list.files(path="./", pattern=".*_family_size_histogram.txt")

### YAML table
### Header
lines <- c(
  "custom_data:",
  "    umi_stats:",
  "        plot_type: 'generalstats'",
  "        pconfig:",
  "            - Mean:",
  "                namespace: 'UMI'",
  "                description: 'Mean number of UMIs per sequence'",
  "                format: '{:,.0f}'",
  "            - Median:",
  "                namespace: 'UMI'",
  "                description: 'Median number of UMIs per sequence'",
  "                format: '{:,.0f}'",
  "            - Standard deviation:",
  "                namespace: 'UMI'",
  "                description: 'Standard deviation of the number of UMIs per sequence'",
  "                format: '{:,.0f}'",
  "            - Min:",
  "                namespace: 'UMI'",
  "                description: 'Minimum number of UMIs per sequence'",
  "                format: '{:,.0f}'",
  "            - Max:",
  "                namespace: 'UMI'",
  "                description: 'Maximum number of UMIs per sequence'",
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
                Mean: %f
                Median: %f
                Standard deviation: %f
                Min: %i
                Max: %i", samp, mean(tostat), median(tostat), sd(tostat), min(tostat), max(tostat))
)
  cat(lines, sep="\n", file="./umi_table_mqc.yaml", append=TRUE)
}
