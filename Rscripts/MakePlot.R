#!/usr/bin/env Rscript

#This is a terminal-executable Rscript to make a plot from the average output across runs

args <- commandArgs(TRUE)
file <- args[1]
outputFile <- args[2]

a <- scan(file, what="character", sep="\n")
b <- sapply(a, strsplit, split=" ")

groups <- NULL
for(i in sequence(length(b))){
  groups <- c(groups, strsplit(b[[i]][1], "\\d+")[[1]][2])
}
groups <- unique(groups)

table <- data.frame(matrix(ncol=4, nrow=length(b)))
colnames(table) <- c("s", "type", "mean", "sd")
for(i in sequence(length(b))){
  table[i,1] <- as.numeric(gsub("[a-z:A-Z]", "", b[[i]][1]))
  #table[i,2] <- groups[agrep(b[[i]][1], groups, max=3)]
  table[i,2] <- strsplit(b[[i]][1], "\\d+")[[1]][2]
  table[i,3] <- as.numeric(gsub("mean=", "", b[[i]][2]))
  table[i,4] <- as.numeric(gsub("sd=", "", b[[i]][3]))
}
table <- table[order(table$s),]

plotmin <- min(table[,3]-table[,4])
plotmax <- max(table[,3]+table[,4])+1

linecols <- rainbow(5)
pointcols <- rainbow(5)

pdf(file=outputFile, width=5, height=5)
plot(table$s, table$mean, type="n", xlab="dataset", ylab="tree length", ylim=c(plotmin,plotmax))
legend("topright", legend=groups, col=linecols, lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
for(i in sequence(length(groups))){
  dataToUse <- table[which(table$type == groups[i]),]
  for(seg in sequence(dim(dataToUse)[1])){
    segments(dataToUse[seg,1], dataToUse[seg,3],dataToUse[seg+1,1], dataToUse[seg+1,3], col=linecols[i])
  }
}
for(i in sequence(length(groups))){
  dataToUse <- table[which(table$type == groups[i]),]
  for(point in sequence(dim(dataToUse)[1])){
    arrows(dataToUse[point,1], dataToUse[point,3]+dataToUse[point,4], dataToUse[point,1], dataToUse[point,3]-dataToUse[point,4], code=3, length=0.05, col=linecols[i], angle=90)
    points(dataToUse[point,1], dataToUse[point,3], bg=pointcols[i], pch=21)
  }
}
dev.off()














