#!/usr/bin/env Rscript

#This is a terminal-executable Rscript to remove invariant sites

library(phrynomics)

args <- commandArgs(TRUE)
file <- args[1]
outputfile <- args[2]
partsFiles <- args[3]

RemoveInvariantSites <- function (SNPdataset, chatty=FALSE, createParts=NULL){
  snpclass <- "table"
  if (class(SNPdataset) == "snp") {
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  snps <- sum(nchar(SNPdataset[1, ]))
  initialLociLengths <- nchar(SNPdataset[1, ])
  splitdata <- SplitSNP(SNPdataset)
  KeepVector <- apply(splitdata, 2, IsVariable)
 
  newSNPdataset <- cSNP(splitdata, KeepVector=KeepVector)
  newsnps <- sum(nchar(newSNPdataset[1, ]))
  if(chatty) 
    message(paste("removed", snps - newsnps, "of", snps, "sites"))
  if(!is.null(createParts)){
    splitInvars <- splitdata[,which(KeepVector == FALSE)]
    invars <- apply(splitInvars, 2, ReturnUniqueBases)
    for(i in sequence(length(invars))){
      if(length(invars[[i]]) > 1)
        invars[[i]] <- invars[[i]][1]  #pull off first invar when there are ambig codes
    }
    invars <- unlist(invars)
    partsInfo <- paste0("ASC_DNA, p1=1-", newsnps)
    invInfo <- paste(length(which(invars == "A")), length(which(invars == "C")), length(which(invars == "T")), length(which(invars == "G")))   #ACTG
    if(createParts == "felsenstein"){
      invInfo <- sum(length(which(invars == "A")), length(which(invars == "C")), length(which(invars == "T")), length(which(invars == "G")))  #ACTG
    }
    if (snpclass == "snp") 
      return(list(data=ReadSNP(newSNPdataset), partsInfo=partsInfo, invInfo=invInfo))
    else
      return(list(data=newSNPdataset, partsInfo=partsInfo, invInfo=invInfo))
  }
  if (snpclass == "snp") 
    return(ReadSNP(newSNPdataset))
  else return(newSNPdataset)
}

RemInvSites <- function(snpFile, outputName, parts=NULL){
  initializeTable <- read.table(snpFile, skip=1, stringsAsFactors=FALSE, row.names=1, colClasses=c("character"))
  snpdata <- ReadSNP(initializeTable)
  snpdata <- RemoveInvariantSites(snpdata, createParts=parts)
  if(!is.null(parts)){
    WriteSNP(snpdata$data, file=outputName)
    pInfo <- paste0("[asc~", outputName, ".inv], ", snpdata$partsInfo)
    write(pInfo, file=paste0(outputName, ".part"))
    write(snpdata$invInfo, file=paste0(outputName, ".inv"))  	
  }
  else
    WriteSNP(snpdata, file=outputName)
}

RemInvSites(file, outputfile, partsFiles)
