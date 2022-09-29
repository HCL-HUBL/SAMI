args <- commandArgs(TRUE)
if(length(args)!=1) stop("USAGE: ./umi_descriptiveStat.R INPUTDIR")
inputdir <- args[1]

library(openxlsx)

### Read all the files
allfile <- list.files(path=inputdir, pattern="_family_size_histogram.txt", full.names=TRUE)

### Get the stats (mean/median/sd/min/max) per sample
res <- NULL
for(ifile in allfile)
{
  tmp    <- read.table(file=ifile, header=TRUE, stringsAsFactors=FALSE)
  samp   <- sub(x=ifile, pattern=".*/", replacement="")
  samp   <- sub(x=samp, pattern="_family_size_histogram.txt", replacement="")
  tostat <- rep(x=tmp$family_size, times=tmp$count)
  res    <- rbind(res,
                  c(Echantillon = samp,
                    Moyenne = mean(tostat),
                    Mediane = median(tostat),
                    "Ecart-type" = sd(tostat),
                    Min = min(tostat),
                    Max = max(tostat)))
}

### Save as CSV and XLSX
write.csv2(x=res, file="./UMI_stats.csv", row.names=FALSE)
write.xlsx(x=res, file="./UMI_stats.xlsx")
