library(phrynomics)

args <- commandArgs(TRUE)
file <- args[1]
lociLength <- as.numeric(args[2])
RElength <- as.numeric(args[3])
KeepRE <- args[4]
outputfile <- args[5]

cSNP51 <- function(splitData, lociLength=51){
#this function will make loci in the lengths of lociLength
  if(class(splitData) == "snp")
    splitData <- splitData$data    
  breaks <- which(1:dim(splitData)[2] %% lociLength == 0)
  breakstart <- 1
  breakend <- breaks[1]
  loci <- data.frame(matrix(nrow=dim(splitData)[1], ncol=length(breaks)))
  rownames(loci) <- rownames(splitData)
  for(i in 1:length(breaks)){
    loci[,i] <- cSNP(splitData[,breakstart:breakend], maintainLoci=FALSE)
    breakstart <- breaks[i]+1
    breakend <- breaks[i+1]
  }
 return(ReadSNP(loci)$data)
}

BreakAtRE <- function(locus, RElength){
#this function takes a locus and just creates a space break at the length of the RE
  splits1 <- cSNP(SplitSNP(as.matrix(locus))[,1:RElength])
  splits2 <- cSNP(SplitSNP(as.matrix(locus))[,(RElength+1):nchar(locus[1])])
  return(cbind(splits1, splits2))
}

GetMostFrequentRE <- function(RE){
  REs <- unique(RE[,1])
  vals <- rep(NA, length(REs))
  for(i in sequence(length(REs))){
    vals[i] <- length(grep(REs[i], RE[,1]))
  }
  return(REs[which.max(vals)])
}

ReplaceWithMissing <- function(RE, MostFreqRE){
#this function will take any infrequent RE and replace with NNNs
  toChange <- which(RE[,1] != MostFreqRE)
  nchars <- nchar(RE[1,])
  RE[toChange,1] <- rep(paste(rep("N", nchars[1]), collapse=""), length(toChange))
  RE[toChange,2] <- rep(paste(rep("N", nchars[2]), collapse=""), length(toChange))
  return(RE)
}

CheckRE <- function(locus, RElength, KeepRE=FALSE){
#this function will take a locus and check if the first X number of sites match (RE), if they do not, then the ones with the lowest frequency are replaced with Ns
  RE <- BreakAtRE(locus, RElength)
  if(length(unique(RE[,1])) == 1){  #if REs are all the same
    if(KeepRE)
      return(locus)  #return original locus (no mods)
    return(RE[,2])  #or return without RE
  }
  if(length(unique(RE[,1])) > 1){
    MostFreqRE <- GetMostFrequentRE(RE)
    newRE <- ReplaceWithMissing(RE, MostFreqRE)
    if(KeepRE)
      return(cSNP(newRE))
    return(newRE[,2])
  }
}

RemNonRandom <- function(snpFile, lociLength, RElength, KeepRE, outputName){
  initializeTable <- read.table(snpFile, skip=1, stringsAsFactors=FALSE, row.names=1, colClasses=c("character"))
  snpdata <- SplitSNP(cSNP(ReadSNP(initializeTable), maintainLoci=FALSE))
  snpsbyloci <- cSNP51(snpdata, lociLength)
  for(i in sequence(dim(snpsbyloci)[2])){
    snpsbyloci[,i] <- CheckRE(snpsbyloci[,i], RElength, KeepRE)
  }
  WriteSNP(snpsbyloci, file=outputName)
}

RemNonRandom(snpFile=file, lociLength, RElength, KeepRE, outputfile)





