#!/usr/bin/env Rscript

#This is a terminal-executable Rscript to 

args <- commandArgs(TRUE)
file <- args[1]
outputFile <- args[2]

a <- scan(file, what="character", sep="\n")
b <- gsub("RAxML_info.", "", a)
c <- gsub(":Tree-Length: ", "_", b)

table <- matrix(ncol=2, nrow=length(c))
for(i in sequence(length(c))){
  table[i,] <- strsplit(c[i], "_")[[1]]
}
uniques <- unique(table[,1])
write("", file=outputFile)
for(i in 1:length(uniques)){
  means <- mean(as.numeric(table[which(table[,1] == uniques[i]), 2]))
  sd <- sd(as.numeric(table[which(table[,1] == uniques[i]), 2]))
  write(paste0(uniques[i], " mean=", means, " sd=", sd), file=outputFile, append=TRUE)
}