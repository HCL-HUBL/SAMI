#!/usr/bin/env Rscript --vanilla

### Read the list of files
allfile <- list.files(path="./", pattern=".*\\.isize\\.txt")

### YAML table
### Header
lines <- c(
  "custom_data:",
  "    insertSize_median:",
  "        plot_type: 'generalstats'",
  "        headers:",
  "            - insertSize.median:",
  "                namespace: 'insertSize'",
  "                description: 'Median insert size'",
  "                format: '{:,.1f}'",
  "        data:"
)
cat(lines, sep="\n", file="./isize_table.yaml")

### Data
for(ifile in allfile)
{
  samp <- gsub(x=ifile, pattern="\\.isize\\.txt", replacement="")
  insertSize <- read.delim(file=ifile, header=FALSE)[,1,drop=TRUE]
  lines <- sprintf("            '%s':
                insertSize.median: %f", samp, median(abs(insertSize)))

  cat(lines, sep="\n", file="./isize_table.yaml", append=TRUE)
}
