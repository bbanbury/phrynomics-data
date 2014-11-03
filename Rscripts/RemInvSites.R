#!/usr/bin/env Rscript

#This is a terminal-executable Rscript to remove invariant sites

library(phrynomics)

args <- commandArgs(TRUE)
file <- args[1]
amount <- args[2]
outputfile <- args[3]


RemoveInvariantSites <- function (SNPdataset, amount=1.0, chatty = FALSE){
  snpclass <- "table"
  if (class(SNPdataset) == "snp") {
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  snps <- sum(nchar(SNPdataset[1, ]))
  initialLociLengths <- nchar(SNPdataset[1, ])
  splitdata <- SplitSNP(SNPdataset)
  KeepVector <- apply(splitdata, 2, IsVariable)
  breaks <- which(splitdata[1, ] == " ")
  if(amount < 1.0){
    invars <- which(KeepVector == "FALSE")
    randomvars <- sample(invars, size=length(invars)*(1-amount))
    KeepVector[randomvars] <- rep("TRUE", length(randomvars))
  }
  newSNPdataset <- cSNP(splitdata, KeepVector = KeepVector, maintainLoci = TRUE)
  newsnps <- sum(nchar(newSNPdataset[1, ]))
  if (chatty) 
    message(paste("removed", snps - newsnps, "of", snps, "sites"))
  if (snpclass == "snp") 
    return(ReadSNP(newSNPdataset))
  else return(newSNPdataset)
}

RemInvSites <- function(snpFile, amount, outputName){
  initializeTable <- read.table(snpFile, skip=1, stringsAsFactors=FALSE, row.names=1, colClasses=c("character"))
  snpdata <- ReadSNP(initializeTable)
  snpdata <- RemoveInvariantSites(snpdata, amount=amount)
  WriteSNP(snpdata, file=outputName)
  #print(snpdata)
}

RemInvSites(file, amount, outputfile)
















