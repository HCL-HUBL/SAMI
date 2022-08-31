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
cat(lines, sep="\n")
