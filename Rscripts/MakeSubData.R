#!/usr/bin/env Rscript

#This is a terminal-executable Rscript create s3,s4,s5 datasets

library(phrynomics)

args <- commandArgs(TRUE)
file <- args[1]
datasets <- args[2]

SNPdataset <- ReadSNP(file)
datasets <- grep(datasets, 1:SNPdataset$ntax):SNPdataset$ntax


coverage <- apply(MakePresentAbsent(SNPdataset$data), 2, sum)

for(i in sequence(length(datasets))){
  keepVec <- coverage >= datasets[i]
  newdata <- as.matrix(SNPdataset$data[,keepVec], nrow=SNPdataset$ntax)
  rownames(newdata) <- rownames(SNPdataset$data)
  WriteSNP(newdata, file=paste0("s", datasets[i], ".txt"))
}