#!/usr/bin/env Rscript --vanilla

##########################
### Get stats on calib ###
##########################

### Read the input data
res <- NULL
for(ifile in list.files(path=".", pattern="_gte[0-9]+\\.cluster$"))
{
	## tmp <- read.table(file=ifile, sep="\t", quote=NULL, comment.char="", colClasses=c("numeric", "numeric", "numeric", NULL, NULL, NULL, NULL, NULL, NULL))
	### Check the file is not empty
	if(file.size(ifile)!=0L)
	{
		tmp <- read.table(file=ifile, sep="\t", quote=NULL, comment.char="")[,1:2]
		tmp <- apply(X=tmp, MARGIN=2, paste0, "_", ifile)
		res <- rbind(res, tmp)
	}
}
colnames(res) <- c("read_cluster_id", "read_node_id")

### Get the number of time UMI are duplicated
res.table <- table(table(res))

### Make a matrix
res.table.mat <- data.frame(res.table)
colnames(res.table.mat) <- c("family_size", "count")
res.table.mat$family_size <- as.numeric(as.character(res.table.mat$family_size))

### Write the result
outname <- sub(x=ifile, pattern="gte[0-9]+\\.cluster$", replacement="family_size_histogram.txt")
write.table(x=res.table.mat, file=outname, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

