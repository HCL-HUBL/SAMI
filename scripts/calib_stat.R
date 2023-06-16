#!/usr/bin/env Rscript --vanilla

##########################
### Get stats on calib ###
##########################

### Define the column names for calib stats
forColName <- c("read_cluster_id",
  "read_node_id")

### Read the input data
res <- NULL
for(ifile in list.files(path=".", pattern=".*_nbr.*\\.cluster"))
{
  ## tmp <- read.table(file=ifile, sep="\t", quote=NULL, comment.char="", colClasses=c("numeric", "numeric", "numeric", NULL, NULL, NULL, NULL, NULL, NULL))
  tmp <- read.table(file=ifile, sep="\t", quote=NULL, comment.char="")[,1:2]
  tmp <- apply(X=tmp, MARGIN=2, paste0, "_", ifile)
  res <- rbind(res, tmp)
}
colnames(res) <- forColName

### Get the number of time UMI are duplicated
res.table <- table(table(res))

### Make a matrix
res.table.mat <- data.frame(res.table)
colnames(res.table.mat) <- c("family_size", "count")
res.table.mat$family_size <- as.numeric(as.character(res.table.mat$family_size))

### Write the result
outname <- gsub(x=ifile, pattern="nbr.*\\.cluster", replacement="family_size_histogram.txt")
write.table(x=res.table.mat, file=outname, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
