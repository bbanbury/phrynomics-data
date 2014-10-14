##  ----------------------------------------  ##
##                                            ##
##               Simulations                  ##
##          edited: 13 Oct 2014               ##
##                                            ##
##  ----------------------------------------  ##



##  ----------------------------------------  ##
##              Data Simulation               ##
##  ----------------------------------------  ##

library(phangorn)
library(phrynomics)
source("phrynomicsFunctions.R")


mainDir <- "/Users/Barb/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper"
ASCDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/ASC.Correction/"
FullDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/FullData/"
RepDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/RepeatRuns/"
SimDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/Simulation/"

setwd(SimDir)

# Unrooted 4-taxon tree 

s <- "(A:0.2, B:0.05, (C:0.2, D:0.05)X:0.05);"
tree <- read.tree(text=s)

#use multilocus dataset that Adam put toghether from 8 nuclear genes off genbank
#rate A <-> C: 1.059151
#rate A <-> G: 2.367842
#rate A <-> T: 0.956396
#rate C <-> G: 0.938012
#rate C <-> T: 3.525699
#rate G <-> T: 1.000000
rates <- as.numeric(c("0.00", "1.059151", "2.367842", "0.956396", "1.059151", "0.00", "0.938012", "3.525699", "2.367842", "0.938012", "0.00", "1.000000", "0.956396", "3.525699", "1.000000", "0.00"))
rateMatrix <- matrix(rates, nrow=4, ncol=4, byrow=TRUE)
rownames(rateMatrix) <- colnames(rateMatrix) <- c("a", "c", "g", "t")

#run simulation
fullsims <- list()
varsims <- list()
for(i in 1:1000){
  simulation <- simSeq(tree, Q=rateMatrix, l=12000, type="DNA", rate=0.10)  #simulate 12000 loci
  fullsims[[i]] <- toupper(as.character(simulation))
  if(length(which(apply(fullsims[[i]], 2, IsVariable))) > 400){  #chop down to 400
    lim <- which(apply(fullsims[[i]], 2, IsVariable))[400]
    fullsims[[i]] <- fullsims[[i]][,1:lim]
    varsims[[i]] <- RemoveInvariantSites(fullsims[[i]])
    if(i %% 10 == 0)
      print(paste(i, "of 1000"))
  }
  else
   stop("too few chars")
}
save(fullsims, varsims, file="sims2.Rdata")


##  ----------------------------------------  ##
##            Simulation Analyses             ##
##  ----------------------------------------  ##


write("", file="invocations.sim.txt")
for(i in sequence(length(fullsims))){
  WriteSNP(fullsims[[i]], file="full.sim", format="phylip")
  run <- paste("raxmlHPC-PTHREADS -T 6 -s full.sim -m GTRCATI -p ", floor(runif(1, min=1, max=10^6)), " -n sim.full", i, sep="")
  write(run, file="invocations.sim.txt", append=TRUE)
  system(command=run)
  if(i %% 10 == 0)
    print(paste(i, "of 1000"))
}
for(i in sequence(length(varsims))){
  WriteSNP(varsims[[i]], file="var.sim", format="phylip")
  run <- paste("raxmlHPC-PTHREADS -T 6 -s var.sim -m ASC_GTRCAT -V -p ", floor(runif(1, min=1, max=10^6)), " -n sim.asc", i, sep="")
  write(run, file="invocations.sim.txt", append=TRUE)
  system(command=run)
  if(i %% 10 == 0)
    print(paste(i, "of 1000"))
}
for(i in sequence(length(varsims))){
  WriteSNP(varsims[[i]], file="var.sim", format="phylip")
  run <- paste("raxmlHPC-PTHREADS -T 6 -s var.sim -m GTRCAT -V -p ", floor(runif(1, min=1, max=10^6)), " -n sim.nonasc", i, sep="")
  write(run, file="invocations.sim.txt", append=TRUE)
  system(command=run)
  if(i %% 10 == 0)
    print(paste(i, "of 1000"))
}



##  ----------------------------------------  ##
##          Missing Data Simulation           ##
##  ----------------------------------------  ##


load("sims2.Rdata")

randomMD <- list()
cdMD <- list()
for(i in sequence(length(varsims))){
  if(i%%10 == 0)
    print(i)
  randomMD[[i]] <- RemoveInvariantSites(deleteData(varsims[[i]], type="random"))
  cdMD[[i]] <- RemoveInvariantSites(deleteData(varsims[[i]], taxa=c("C", "D")))
}
save(randomMD, cdMD, file="missingData.Rdata")


##  ----------------------------------------  ##
##           Missing Data Analyses            ##
##  ----------------------------------------  ##

write("", file="invocations.md.txt")
for(i in sequence(length(randomMD))){
  WriteSNP(randomMD[[i]], file="randMD.sim", format="phylip")
  run <- paste("raxmlHPC-PTHREADS -T 6 -s randMD.sim -m ASC_GTRCAT -V -p ", floor(runif(1, min=1, max=10^6)), " -n randMD.asc", i, sep="")
  write(run, file="invocations.md.txt", append=TRUE)
  system(command=run)
  run <- paste("raxmlHPC-PTHREADS -T 6 -s randMD.sim -m GTRCAT -V -p ", floor(runif(1, min=1, max=10^6)), " -n randMD.nonasc", i, sep="")
  write(run, file="invocations.md.txt", append=TRUE)
  system(command=run)
  if(i %% 10 == 0)
    print(paste(i, "of 1000"))
}
for(i in sequence(length(cdMD))){
  WriteSNP(cdMD[[i]], file="cdMD.sim", format="phylip")
  run <- paste("raxmlHPC-PTHREADS -T 6 -s cdMD.sim -m ASC_GTRCAT -V -p ", floor(runif(1, min=1, max=10^6)), " -n cdMD.asc", i, sep="")
  write(run, file="invocations.md.txt", append=TRUE)
  system(command=run)
  run <- paste("raxmlHPC-PTHREADS -T 6 -s cdMD.sim -m GTRCAT -V -p ", floor(runif(1, min=1, max=10^6)), " -n cdMD.nonasc", i, sep="")
  write(run, file="invocations.md.txt", append=TRUE)
  system(command=run)
  if(i %% 10 == 0)
    print(paste(i, "of 1000"))
}

setwd(mainDir)




















