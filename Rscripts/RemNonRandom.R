library(phrynomics)

args <- commandArgs(TRUE)
file <- args[1]
rootSeqfile <- args[2]
lociLength <- as.numeric(args[3])
RElength <- as.numeric(args[4])
KeepRE <- args[5]
outputfile <- args[6]

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

#GetMostFrequentRE <- function(RE){
#  REs <- unique(RE[,1])
#  vals <- rep(NA, length(REs))
#  for(i in sequence(length(REs))){
#    vals[i] <- length(grep(REs[i], RE[,1]))
#  }
#  return(REs[which.max(vals)])
#}

GetShortRE <- function(rootlocus, RElength){
  Rep <- strsplit(gsub(" ", "", rootlocus), "")[[1]][1: RElength]
  return(paste(Rep, collapse=""))
}

GetREsFromRoot <- function(rootSeqfile){
  lines <- system(paste("grep 'root seq'", rootSeqfile), intern=TRUE)
  lines <- gsub("locus  51 root seq              ", "", lines)
  REs <- sapply(lines, GetShortRE, RElength=12, USE.NAMES=FALSE)
  return(REs)
}

ReplaceWithMissing <- function(RE, rootRE){
#this function will take any infrequent RE and replace with NNNs
  toChange <- which(RE[,1] != rootRE)
  nchars <- nchar(RE[1,])
  RE[toChange,1] <- rep(paste(rep("N", nchars[1]), collapse=""), length(toChange))
  RE[toChange,2] <- rep(paste(rep("N", nchars[2]), collapse=""), length(toChange))
  return(RE)
}

CheckRE <- function(locus, rootRE, RElength, KeepRE=FALSE){
#This function will check that a locus matches a RE, if it does not, it replaces with NNNs
  RE <- BreakAtRE(locus, RElength)
  if(any(rootRE != RE[,1])){
    dontmatch <- which(rootRE != RE[,1])
    RE <- ReplaceWithMissing(RE, rootRE)
  }
  if(KeepRE)
    return(cSNP(RE))
  return(RE[,2])
}

RemNonRandom <- function(snpFile, rootSeqfile, lociLength, RElength, KeepRE, outputName){
  initializeTable <- read.table(snpFile, skip=1, stringsAsFactors=FALSE, row.names=1, colClasses=c("character"))
  snpdata <- SplitSNP(cSNP(ReadSNP(initializeTable), maintainLoci=FALSE))
  snpsbyloci <- cSNP51(snpdata, lociLength)
  rootREs <- GetREsFromRoot(rootSeqfile)
  for(i in sequence(dim(snpsbyloci)[2])){
    snpsbyloci[,i] <- CheckRE(snpsbyloci[,i], rootREs[i], RElength, KeepRE)
  }
  WriteSNP(snpsbyloci, file=outputName)
}

RemNonRandom(snpFile=file, rootSeqfile, lociLength, RElength, KeepRE, outputfile)





