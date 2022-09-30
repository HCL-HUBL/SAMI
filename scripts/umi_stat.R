args <- commandArgs(TRUE)
if(length(args)!=2) stop("USAGE: ./umi_stat.R SAMPLE_ID INPUT_FILE")
sample    <- args[1]
inputFile <- args[2]

if(!file.exists(inputFile)) stop("INPUT_FILE must exist")

### Read the file
umi_hist <- read.delim(file=inputFile)

# YAML graph
lines <- c(
  "id: 'UMI_duplication'",
  "section_name: 'UMI duplication'",
  "description: 'distribution of the number of UMIs depending on their number of copies.'",
  "plot_type: 'linegraph'",
  "pconfig:",
  "    id: 'UMI_duplication_linegraph'",
  "    title: 'UMI duplication'",
  "    xlab: 'Number of UMI copies'",
  "    ylab: 'Fraction of the total number of read'",
  "data:",
  sprintf("    %s: { %s }", sample, paste(sprintf("%i: %g", umi_hist$family_size, umi_hist$fraction), collapse=", "))
)
cat(lines, sep="\n", file=sprintf("./%s_mqc.yaml", sample))

# YAML table
tostat <- rep(x=umi_hist$family_size, times=umi_hist$count)
lines <- c(
  "custom_data:",
  "\tumi_stats:",
  "\t\tplot_type: 'generalstats'",
  "\t\tpconfig:",
  "\t\t\t- Mean:",
  "\t\t\t\tnamespace: 'UMI'",
  "\t\t\t\tdescription: 'Mean number of UMIs per sequence'",
  "\t\t\t\tformat: '{:,.0f}'",
  "\t\t\t- Median:",
  "\t\t\t\tnamespace: 'UMI'",
  "\t\t\t\tdescription: 'Median number of UMIs per sequence'",
  "\t\t\t\tformat: '{:,.0f}'",
  "\t\t\t- Standard deviation:",
  "\t\t\t\tnamespace: 'UMI'",
  "\t\t\t\tdescription: 'Standard deviation of the number of UMIs per sequence'",
  "\t\t\t\tformat: '{:,.0f}'",
  "\t\t\t- Min:",
  "\t\t\t\tnamespace: 'UMI'",
  "\t\t\t\tdescription: 'Minimum number of UMIs per sequence'",
  "\t\t\t\tformat: '{:,.0f}'",
  "\t\t\t- Max:",
  "\t\t\t\tnamespace: 'UMI'",
  "\t\t\t\tdescription: 'Maximum number of UMIs per sequence'",
  "\t\t\t\tformat: '{:,.0f}'",
  "\t\tdata:",
  sprintf("\t\t\t'%s':
      \t\t\tMean: %f
      \t\t\tMedian: %f
      \t\t\tStandard deviation: %f
      \t\t\tMin: %i
      \t\t\tMax: %i", sample, mean(tostat), median(tostat), sd(tostat), min(tostat), max(tostat))
)
cat(lines, sep="\n", file=sprintf("./%s_table_mqc.yaml", sample))

## tostat <- rep(x=umi_hist$family_size, times=umi_hist$count)
## lines <- c(
##   "id: 'UMI_stat'",
##   "section_name: 'UMI stat'",
##   "description: 'descriptive statistics on the number of UMIs.'",
##   "plot_type: 'generalstats'",
##   "pconfig:",
##   "    id: 'UMI_stat_generalstats'",
##   "    title: 'UMI stat'",
##   "data:",
##   sprintf("%s: {
##       Moyenne : %s,
##       Mediane : %s,
##       Ecart-type : %s,
##       Min : %s,
##       Max : %s
## }", sample, mean(tostat), median(tostat), sd(tostat), min(tostat), max(tostat))
## )
## cat(lines, sep="\n", file=sprintf("./%s_table_mqc.yaml", sample))
