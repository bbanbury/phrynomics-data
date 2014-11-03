#!/usr/bin/env Rscript

#This is a terminal-executable Rscript to remove mising data

library(phrynomics)

args <- commandArgs(TRUE)
file <- args[1]
amount <- as.numeric(args[2])
tax <- args[3]
outputfile <- args[4]

cSNP39 <- function(splitData, lociLength=39){
#this function will make 39 base loci
  if(class(splitData) == "snp")
    splitData <- splitData$data
  breaks <- which(1:dim(splitData)[2] %% lociLength == 0)
  breakstart <- 1
  breakend <- breaks[1]
  loci <- data.frame(matrix(nrow=dim(splitData)[1], ncol=length(breaks)))
  rownames(loci) <- rownames(splitData)
  for(i in 1:length(breaks)){
    loci[,i] <- cSNP(splitData[,breakstart:breakend])
    breakstart <- breaks[i]+1
    breakend <- breaks[i+1]
  }
 return(ReadSNP(loci)$data)
}

deleteData <- function(data, taxa=NULL, percent=0.4){
#function to remove loci from certain taxa
  if(class(data) == "snp")
    data <- data$data
  whichTax <- grep(taxa, rownames(data))
  data <- as.matrix(SplitSNP(cSNP(data)))
  data <- cSNP39(data)
  chars <- paste(rep.int("N", 39), collapse="")
  toDel <- sample(dim(data)[2], dim(data)[2]*percent)
  data[whichTax, toDel] <- rep(chars, length(toDel)*length(whichTax))
  return(as.data.frame(data, stringsAsFactors=FALSE))
}

RemMissing <- function(snpFile, snpAmount, snpTax, outputName){
  initializeTable <- read.table(snpFile, skip=1, stringsAsFactors=FALSE, row.names=1, colClasses=c("character"))
  snpdata <- ReadSNP(initializeTable)
  remdata <- deleteData(data=snpdata, taxa=snpTax, percent=snpAmount)
  WriteSNP(remdata, file=outputName)
}

RemMissing(snpFile=file, amount, tax, outputfile)








