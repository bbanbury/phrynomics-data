##  ----------------------------------------  ##
##                                            ##
##           Phrynomics Analyses              ##
##           edited: 27 Oct 2014              ##
##                                            ##
##  ----------------------------------------  ##

library(phangorn)
library(phrynomics)
source("~/phrynomics-data/trunk/phrynomicsFunctions.R")

mainDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/cResults-k80"
FigDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/TablesFigures-k80"



##  ----------------------------------------  ##
##       Post-Analyses Trees and Data         ##
##  ----------------------------------------  ##

#  Load RAxML Trees and post-analyses scraping

setwd(mainDir)
fullfiles <- system("ls s*full.phy", intern=T)
snpfiles <- system("ls s*noAmbigs.phy", intern=T)


#load full, vstam, vfel, vlewis
analysis <- "RAxML"
RAxML.trees <- system("ls RAxML_bipartitions.*", intern=T)  #trees with bootstraps
RAxML.TreeList <- CreateTreeList(RAxML.trees, "RAxML")
RAxML.TreeList <- lapply(RAxML.TreeList, ladderize)
ML.results <- GetRAxMLStatsPostAnalysis(".")

treeMatrix <- CreateTreeMatrix(RAxML.trees)

# full - nonASC
whichRuns <- c("full", "nonasc")
full.nonASC.treeMatrix <- AddBLD(AddTreeDist(treeMatrix[which(colnames(treeMatrix) %in% whichRuns)], RAxML.TreeList), RAxML.TreeList)
fullnonASC.branchComp <- list()
for(i in sequence(dim(full.nonASC.treeMatrix)[1])) {
  tree1 <- assTrees(full.nonASC.treeMatrix[i, which(colnames(full.nonASC.treeMatrix) == whichRuns[1])], RAxML.TreeList)[[1]]
  tree2 <- assTrees(full.nonASC.treeMatrix[i, which(colnames(full.nonASC.treeMatrix) == whichRuns[2])], RAxML.TreeList)[[1]]
  fullnonASC.branchComp[[i]] <- MakeBranchLengthMatrix(tree1, tree2, dataset=rownames(full.nonASC.treeMatrix)[i])
  names(fullnonASC.branchComp)[[i]] <- rownames(full.nonASC.treeMatrix)[i]
}

# full - lewis
whichRuns <- c("full", "lewis")
full.lewis.treeMatrix <- AddBLD(AddTreeDist(treeMatrix[which(colnames(treeMatrix) %in% whichRuns)], RAxML.TreeList), RAxML.TreeList)
fulllewis.branchComp <- list()
for(i in sequence(dim(full.lewis.treeMatrix)[1])) {
  tree1 <- assTrees(full.lewis.treeMatrix[i, which(colnames(full.lewis.treeMatrix) == whichRuns[1])], RAxML.TreeList)[[1]]
  tree2 <- assTrees(full.lewis.treeMatrix[i, which(colnames(full.lewis.treeMatrix) == whichRuns[2])], RAxML.TreeList)[[1]]
  fulllewis.branchComp[[i]] <- MakeBranchLengthMatrix(tree1, tree2, dataset=rownames(full.lewis.treeMatrix)[i])
  names(fulllewis.branchComp)[[i]] <- rownames(full.lewis.treeMatrix)[i]
}

# full - stam
whichRuns <- c("full", "stam")
full.stam.treeMatrix <- AddBLD(AddTreeDist(treeMatrix[which(colnames(treeMatrix) %in% whichRuns)], RAxML.TreeList), RAxML.TreeList)
fullstam.branchComp <- list()
for(i in sequence(dim(full.stam.treeMatrix)[1])) {
  tree1 <- assTrees(full.stam.treeMatrix[i, which(colnames(full.stam.treeMatrix) == whichRuns[1])], RAxML.TreeList)[[1]]
  tree2 <- assTrees(full.stam.treeMatrix[i, which(colnames(full.stam.treeMatrix) == whichRuns[2])], RAxML.TreeList)[[1]]
  fullstam.branchComp[[i]] <- MakeBranchLengthMatrix(tree1, tree2, dataset=rownames(full.stam.treeMatrix)[i])
  names(fullstam.branchComp)[[i]] <- rownames(full.stam.treeMatrix)[i]
}




#  Create global objects for tables and figs to be made

levels <-NULL
for(i in sequence(length(snpfiles))){
  levels <- as.numeric(c(levels, strsplit(snpfiles[i], "\\D+")[[1]][2])) #numerical datasets from file names
}
orderedLevels <- sort(levels)  #numerical datasets from file names
whichDatasets <- paste("s", levels, "", sep="")  #datasets from file names
AllOrder <- paste("s", seq(5, 70, 5), "", sep="") #datasets from sequence
orderToGo <- AllOrder[AllOrder %in% whichDatasets]  #datasets of sequence that exist
focalDatasets <- c("s5", "s25", "s50")
treeFigDatasets <- c("s5", "s15", "s25", "s35", "s45", "s55", "s65")

save(fullfiles, snpfiles, RAxML.trees, RAxML.TreeList, treeMatrix, ML.results, fullnonASC.branchComp, fulllewis.branchComp, fullstam.branchComp, orderedLevels, orderToGo, focalDatasets, treeFigDatasets, file="phrynoResults.Rdata")



##  ----------------------------------------  ##
##      End Post-Analyses Trees and Data      ##
##  ----------------------------------------  ##






##  ----------------------------------------  ##
##                 Make Tables                ##
##  ----------------------------------------  ##

#  If you want to start from here, you can load up the workspace with all the data.
load("phrynoResults.Rdata")

#  Table 1 was made by hand


#  Make Table 2. Summary ddRadSeq data. 

toGo <- sort(orderedLevels, decreasing=TRUE)
snpdatafiles <- system("ls s*dataonly.snps", intern=TRUE)  #analyses run with concatenated file missing loci info
full.ML.results <- ML.results[which(ML.results[,2] == "full"),]
table2 <- matrix(nrow=length(toGo), ncol=6)
rownames(table2) <- paste("s", toGo, sep="")
colnames(table2) <- c("Matrix", "Loci", "VariableSites", "Missing(%)", "PhrynoOverlap", "non-PhrynoOverlap")
table2[,1] <- paste("s", toGo, sep="")
rows <- match(full.ML.results[,1], rownames(table2))
table2 <- table2[!is.na(rows),]  
table2[,c(2,3,4)] <- as.matrix(full.ML.results[order(rows[which(!is.na(rows))]), c(3,4,6)])
dataOverlap <- list()
for(i in sequence(length(snpdatafiles))){
  dataset <- ReadSNP(snpdatafiles[[i]], fileFormat="phy", extralinestoskip=1)
  datasetName <- strsplit(strsplit(snpdatafiles[i], "/")[[1]][length(strsplit(snpdatafiles[i], "/")[[1]])], "dataonly")[[1]][1]
  table2[which(datasetName == rownames(table2)), 2] <- dataset$nloci
  dataOverlap[[i]] <- DataOverlap(dataset)[[2]]
  names(dataOverlap)[[i]] <- datasetName
}
for(m in sequence(dim(table2)[1])){
  whichDataset <- which(names(dataOverlap) == rownames(table2)[m])
  table2[m,6] <- round(mean(dataOverlap[[whichDataset]][-grep(pattern="PH", names(dataOverlap[[whichDataset]]))]), digits=2)
  table2[m,5] <- round(mean(dataOverlap[[whichDataset]][grep(pattern="PH", names(dataOverlap[[whichDataset]]))]), digits=2)
}
setwd(FigDir)
write.table(table2, file="table2.txt", quote=FALSE, row.names=FALSE)  


#  Table 3. Support for Clades

subtrees <- RAxML.TreeList[grep("s50|s25|s5", names(RAxML.TreeList))]
subtrees <- subtrees[-grep("s55", names(subtrees))]  #grep pulls s50 with s5
subtreeLevels <- paste0("s", c(50, 25, 5))
trees50 <- subtrees[grep("s50", names(subtrees))]
trees25 <- subtrees[grep("s25", names(subtrees))]
trees5 <- subtrees[grep("s5", names(subtrees))]
trees5 <- trees5[-grep("s50", names(trees5))]

Sceloporus <- c("SCMA1", "SCOC1", "SCGA1", "SCAN1")
Sceloporinae <- c(Sceloporus, "UROR1", "URBI1", "UTST1", "PETH1")
Tapaja <- c("PHHE5", "PHHE4", "PHHE2", "PHHE3", "PHHE1", "PHDI2", "PHDI1", "PHDO1", "PHDO2", "PHOR1", "PHOR2", "PHOR3", "PHOR4")
Doliosaurus <- c("PHMO1", "PHMO2", "PHGO3", "PHGO2", "PHGO1", "PHGO4", "PHPL1", "PHPL2", "PHPL3")
Brevicauda <- c("PHSH1", "PHSH4", "PHSH2", "PHSH3", "PHTA2", "PHTA3", "PHTA1", "PHTA4", "PHBR1", "PHBR3", "PHBR2", "PHBR4")
Anota <- c("PHMC1", "PHMC3", "PHMC2", "PHMC4", "PHSO3", "PHSO2", "PHSO1", "PHCE1", "PHCE3", "PHCE2", "PHCE4", "PHBL3", "PHBL2", "PHBL1", "PHBL4", "PHCE6", "PHCO1", "PHCE5")
Phrynosomatini <- c("PHAS4", "PHAS1", "PHAS3", "PHAS2", Anota, "PHCN1", "PHCN2", "PHCN3", "PHCN4", Brevicauda, Doliosaurus, Tapaja)
Callisaurini <- c("COTE1", "HOMA1", "CADR1", "CADR2", "UMNO1")
Phrynosomatinae <- c(Callisaurini, Phrynosomatini)
HolbrookiaCalli <- c("HOMA1", "CADR1", "CADR2")
HolbrookiaCopho <- c("HOMA1", "COTE1")

clades <- list(Sceloporinae=Sceloporinae, Phrynosomatinae=Phrynosomatinae, Callisaurini=Callisaurini, HolbrookiaCalli=HolbrookiaCalli, HolbrookiaCopho=HolbrookiaCopho, Phrynosomatini=Phrynosomatini, Anota=Anota, Brevicauda=Brevicauda, Doliosaurus=Doliosaurus, Tapaja=Tapaja, Sceloporus=Sceloporus)

table3 <- matrix(nrow=length(clades), ncol=4)
rownames(table3) <- names(clades)
colnames(table3) <- c("full", "nonasc", "lewis", "stam")

for(i in sequence(dim(table3)[1])){
  taxa <- clades[[which(names(clades) == rownames(table3)[i])]]
  for(col in sequence(dim(table3)[2])){
    colnam <- colnames(table3)[col]
    threeVals <- rep("--", 3)
    tree50 <- trees50[[grep(colnam, names(trees50))]]
    tree25 <- trees25[[grep(colnam, names(trees25))]]
    tree5 <- trees5[[grep(colnam, names(trees5))]]
    trees <- c(tree50, tree25, tree5)
    for(num in 1:3){
      tree <- trees[[num]]
      mrca <- getMRCA(tree, taxa)
      t2 <- nodeLeaves(tree, mrca)
      if(all(taxa %in% t2) && all(t2 %in% taxa)){
        EL <- GetEdgeList(tree)
        sup <- EL[which(EL[,2] == mrca), 5]
        threeVals[num] <- sup
      }
    }
    table3[i,col] <- paste(threeVals, collapse="/")
  }
}
write.table(table3, file="table3.txt", quote=FALSE, sep=" & ")
 


# Table 4. RF distances. 

RFs <- cbind(full.nonASC.treeMatrix$symmetric.difference, full.lewis.treeMatrix$symmetric.difference, full.stam.treeMatrix$symmetric.difference )
rownames(RFs) <- rownames(full.lewis.treeMatrix)
colnames(RFs) <- c("nonasc", "lewis", "stam") 
write.table(RFs, file="RFs.txt")


# RF dists between full comparisons
RFdist <- rep(NA, 13)
fullTrees <- RAxML.TreeList[grep("full", names(RAxML.TreeList))]
for(i in sequence(length(orderToGo)-1)){
  tree1 <- fullTrees[[grep(paste0(orderToGo[i], "full"), names(fullTrees))]]
  tree2 <- fullTrees[[grep(paste0(orderToGo[i+1], "full"), names(fullTrees))]]
  RFdist[i] <- phangorn::treedist(tree1, tree2)[[1]]
  names(RFdist)[i] <- paste0(orderToGo[i], "-", orderToGo[i+1])
}
write(RFdist, file="fullRFdists.txt")


##  ----------------------------------------  ##
##               End Make Tables              ##
##  ----------------------------------------  ##






##  ----------------------------------------  ##
##                Make Figures                ##
##  ----------------------------------------  ##


#  Figure 1 -- Simulation, done using Rscripts


#  Figure 2 -- Phylogeny for group, made by hand


#  Make Figure 3. Tree length comparisons

setwd(FigDir)
pdf(file="Figure2.pdf", width=5, height=5)
cols <- rainbow(length(unique(ML.results$Model)))
plot(rep(orderedLevels, length(unique(ML.results$Model))), ML.results$TreeLength, type="n", ylab="Tree Length", xlab="Missing Data")
title(main="RAxML")
legend("topright", legend=unique(ML.results$Model), col=cols, lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
for(i in sequence(length(orderToGo)-1)){
  dataToUse <- which(orderToGo[i] == ML.results$Level)
  nextDataToUse <- which(orderToGo[i+1] == ML.results$Level)
  for(j in sequence(unique(ML.results$Model))){
    segments(orderedLevels[i], ML.results$TreeLength[dataToUse[j]], orderedLevels[i+1], ML.results$TreeLength[nextDataToUse[j]], col=cols[j])
  }
}
for(i in sequence(length(orderToGo))){
  dataToUse <- which(orderToGo[i] == ML.results$Level)
  for(j in sequence(length(unique(ML.results$Model)))){
    points(orderedLevels[i], ML.results$TreeLength[dataToUse[j]], pch=21, bg="white")
  }
}
dev.off()


#  Make Figure 4. Scatterplot branch lengths

setwd(FigDir)
usr2dev <- function(x) 180/pi*atan(x*diff(par("usr")[1:2])/diff(par("usr")[3:4]))
getX <- function(y, linmod) (y-linmod$coefficients[1])/linmod$coefficients[2]  #(y -b)/m = x
pdf(file="Figure3.pdf", width=8.5, height=5)
op <- par(mar=par("mar")/1.7)
layout(matrix(1:6, nrow=2, byrow=TRUE), respect=TRUE)  
analyses <- focalDatasets[c(3,1)]
for(anal in 1:length(analyses)){
  whichAnalysis <- analyses[anal] 
  data1 <- which(whichAnalysis == names(fullnonASC.branchComp))
  data2 <- which(whichAnalysis == names(fulllewis.branchComp))
  data3 <- which(whichAnalysis == names(fullstam.branchComp))
  BLs1 <- fullnonASC.branchComp[[data1]]$branchlength[which(fullnonASC.branchComp[[data1]]$present)]
  corr.BLs1 <- fullnonASC.branchComp[[data1]]$corr.BL[which(fullnonASC.branchComp[[data1]]$present)]
  BLs2 <- fulllewis.branchComp[[data2]]$branchlength[which(fulllewis.branchComp[[data2]]$present)]
  corr.BLs2 <- fulllewis.branchComp[[data2]]$corr.BL[which(fulllewis.branchComp[[data2]]$present)]
  BLs3 <- fullstam.branchComp[[data3]]$branchlength[which(fullstam.branchComp[[data3]]$present)]
  corr.BLs3 <- fullstam.branchComp[[data3]]$corr.BL[which(fullstam.branchComp[[data3]]$present)]
  BLs <- list(nonASC=BLs1, lewis=BLs2, stam=BLs3)    
  corr.BLs <- list(nonASC=corr.BLs1, lewis=corr.BLs2, stam=corr.BLs3)
  for(plotty in sequence(length(BLs))){
    #print(paste("you are at", plotty))
    maxlims <- max(BLs[[plotty]], corr.BLs[[plotty]])
    plot(BLs[[plotty]], corr.BLs[[plotty]], pch=21, bg="gray", xlab="full", ylab=names(BLs)[plotty], xlim=c(0, maxlims), ylim=c(0, maxlims), type="n")
    linmod <- lm(corr.BLs[[plotty]] ~ BLs[[plotty]])
    abline(linmod, lty=2)
    y <- 0.18
    points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
    points(BLs[[plotty]], corr.BLs[[plotty]], pch=21, bg="gray")
    text(x=getX(y, linmod), y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
    lines(c(-1,1), c(-1,1))
  }
  title(main=analyses[anal])
}
dev.off()


#  Figure 5a. Colored branch SNP phylogenies -- All vs lewis

setwd(FigDir)
pdf(file="Figure4a.pdf", width=11, height=5)
par(mar = c(1,0,1,0))
layout(matrix(1:7, nrow=1, byrow=TRUE), respect=TRUE)
for(fig in treeFigDatasets){
  whichAnalysis <- fig 
  print(whichAnalysis)
  dataToUse <- which(rownames(treeMatrix) == fig)
  BL.AllTrees <- fulllewis.branchComp
  tree1 <- assTrees(treeMatrix[dataToUse, grep("full", colnames(treeMatrix))], RAxML.TreeList)[[1]]
  tree2 <- assTrees(treeMatrix[dataToUse, grep("lewis", colnames(treeMatrix))], RAxML.TreeList)[[1]]
  legtxt <- c("> 500%", "> 400%", "> 300%", "> 200%", "> 100%", "\u00b1 100%", "< -100%")
  legcolors <- c(rgb(153,0,0, max=255), rgb(225,0,0, max=255), rgb(225,128,0, max=255), rgb(255,178,102, max=255), rgb(255,255,102, max=255), "gray", rgb(51,51,255, max=255))
  edgeColors <- BL.AllTrees[[dataToUse]]$edgeColor1
  if(length(edgeColors[which(is.na(edgeColors))]) > 0)
    edgeColors[which(is.na(edgeColors))] <- rep("gray", length(which(is.na(edgeColors))))
  print(paste(names(BL.AllTrees)[[dataToUse]], mean(BL.AllTrees[[dataToUse]]$relativeBLdiff, na.rm=TRUE)))
  plot(tree1, edge.lty=BL.AllTrees[[dataToUse]]$edgelty, edge.color= edgeColors, cex=0.5, edge.width=1, show.tip.label=FALSE)
  if(fig == "s5")
    legend("bottomleft", legend=legtxt, col=legcolors, lty=rep(1,7), lwd=2, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Branchlength Difference")))
  title(main= whichAnalysis)
}
dev.off()


#  Figure 5b. Colored branch SNP phylogenies -- All vs stam

setwd(FigDir)
pdf(file="Figure4b.pdf", width=11, height=5)
par(mar = c(1,0,1,0))
layout(matrix(1:7, nrow=1, byrow=TRUE), respect=TRUE)
for(fig in treeFigDatasets){
  whichAnalysis <- fig 
  #print(whichAnalysis)
  dataToUse <- which(rownames(treeMatrix) == fig)
  BL.AllTrees <- fullstam.branchComp
  tree1 <- assTrees(treeMatrix[dataToUse, grep("full", colnames(treeMatrix))], RAxML.TreeList)[[1]]
  tree2 <- assTrees(treeMatrix[dataToUse, grep("stam", colnames(treeMatrix))], RAxML.TreeList)[[1]]
  legtxt <- c("> 25%", "\u00b1 25%", "< 25%", "< 50%", "< 100%")
  legcolors <- c(rgb(196,156,100, max=255), "gray", rgb(51,255,241, max=255), rgb(51,194,255, max=255), rgb(51,51,255, max=255))  
  edgeColors <- BL.AllTrees[[dataToUse]]$edgeColor2
  if(length(edgeColors[which(is.na(edgeColors))]) > 0)
    edgeColors[which(is.na(edgeColors))] <- rep("gray", length(which(is.na(edgeColors))))
  #print(paste(names(BL.AllTrees)[[dataToUse]], mean(BL.AllTrees[[dataToUse]]$relativeBLdiff, na.rm=TRUE)))
  plot(tree1, edge.lty=BL.AllTrees[[dataToUse]]$edgelty, edge.color= edgeColors, cex=0.5, edge.width=1, show.tip.label=FALSE)
  if(fig == "s5")
    legend("bottomleft", legend=legtxt, col=legcolors, lty=rep(1,7), lwd=1, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Branchlength Difference")))
  title(main= whichAnalysis)
}
dev.off()



#  Make Figure 6. Scatterplot branch length support 

setwd(FigDir)
pdf(file="Figure6a-d.pdf", width=8.5, height=5)
par(mar=rep(4,4))
layout(matrix(1:6, nrow=2, byrow=TRUE), respect=TRUE)  
analyses <- focalDatasets[c(3,1)]
for(whichAnalysis in focalDatasets[c(3,1)]){
  data1 <- which(whichAnalysis == names(fullnonASC.branchComp))
  data2 <- which(whichAnalysis == names(fulllewis.branchComp))
  data3 <- which(whichAnalysis == names(fullstam.branchComp))
  data1rows <- which(fullnonASC.branchComp[[data1]]$present)[which(fullnonASC.branchComp[[data1]]$present) %in% which(fullnonASC.branchComp[[data1]]$support != 0)]  #don't want all the tip support 0s or non homologous clades
  data2rows <- which(fulllewis.branchComp[[data2]]$present)[which(fulllewis.branchComp[[data2]]$present) %in% which(fulllewis.branchComp[[data2]]$support != 0)]  #don't want all the tip support 0s or non homologous clades
  data3rows <- which(fullstam.branchComp[[data3]]$present)[which(fullstam.branchComp[[data3]]$present) %in% which(fullstam.branchComp[[data3]]$support != 0)]  #don't want all the tip support 0s or non homologous clades
  support1 <- fullnonASC.branchComp[[data1]]$support[data1rows]
  corr.support1 <- fullnonASC.branchComp[[data1]]$corr.support[data1rows]
  support2 <- fulllewis.branchComp[[data2]]$support[data2rows]
  corr.support2 <- fulllewis.branchComp[[data2]]$corr.support[data2rows]
  support3 <- fullstam.branchComp[[data3]]$support[data3rows]
  corr.support3 <- fullstam.branchComp[[data3]]$corr.support[data3rows]
  xlims <- ylims <- c(1,100)
  supports <- list(nonasc=support1, lewis=support2, stam=support3)
  corr.supports <- list(nonasc=corr.support1, lewis=corr.support2, stam=corr.support3)
  for(plotty in sequence(length(supports))){
    plot(supports[[plotty]], corr.supports[[plotty]], ylab=names(corr.supports)[plotty], xlab="full", pch=21, bg="gray", xlim=xlims, ylim=ylims)  
    title(main=whichAnalysis)
    linmod <- lm(corr.supports[[plotty]] ~ supports[[plotty]])
    abline(linmod, lty=2)
    text(x=80, y=10, paste("r =", round(cor(corr.supports[[plotty]], supports[[plotty]]), digits=2)))
  }
}
dev.off()


#  Figure 7a. -- Time

setwd(mainDir)
infoFiles <- system("ls RAxML_info*", intern=TRUE)
infoFiles <- infoFiles[-grep("nonasc", infoFiles)]
timeMatrix <- CreateTimeMatrix(infoFiles)
setwd(FigDir)
pdf(file="Figure7.pdf", width=5, height=5)
cols <- rainbow(dim(timeMatrix)[2])
plot(rep(orderedLevels, dim(timeMatrix)[2]), unlist(timeMatrix), type="n", ylab="Hours", xlab="Dataset", axes=FALSE)
axis(side=2)
axis(side=1, at=orderedLevels, labels=paste0("s", orderedLevels))
box()
index <- 0
title(main="Time")
for(col in colnames(timeMatrix)){
  index <- index+1
  for(i in sequence(length(orderToGo)-1)){
    dataToUse <- which(orderToGo[i] == rownames(timeMatrix))
    nextDataToUse <- which(orderToGo[i+1] == rownames(timeMatrix))
    segments(orderedLevels[i], timeMatrix[dataToUse, col], orderedLevels[i+1], timeMatrix[nextDataToUse, col], col=cols[index])
  }
}
legend("topright", legend=colnames(timeMatrix), col=cols, lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
dev.off()


#  Figure 7b. -- Site Patterns
setwd(mainDir)
infoFiles <- system("ls RAxML_info*", intern=TRUE)
spMatrix <- CreateSitePatternMatrix(infoFiles)
setwd(FigDir)
pdf(file="Figure7b.pdf", width=5, height=5)
cols <- rainbow(dim(spMatrix)[2])
plot(rep(orderedLevels, dim(spMatrix)[2]), unlist(spMatrix), type="n", ylab="Site Patterns", xlab="Dataset", axes=FALSE)
axis(side=2)
axis(side=1, at=orderedLevels, labels=paste0("s", orderedLevels))
box()
index <- 0
title(main="")
for(col in colnames(spMatrix)){
  index <- index+1
  for(i in sequence(length(orderToGo)-1)){
    dataToUse <- which(orderToGo[i] == rownames(spMatrix))
    nextDataToUse <- which(orderToGo[i+1] == rownames(spMatrix))
    segments(orderedLevels[i], spMatrix[dataToUse, col], orderedLevels[i+1], spMatrix[nextDataToUse, col], col=cols[index])
  }
}
legend("topright", legend=colnames(spMatrix), col=cols, lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
dev.off()
write.table(spMatrix, file="SitePatterns.txt")



#  Figure 8. -- RF distances.

setwd(FigDir)
pdf(file="Figure8.pdf", width=5, height=5)
cols <- rainbow(dim(RFs)[2])
index <- 0
plot(rep(levels, dim(RFs)[2]), unlist(RFs), type="n")
for(col in colnames(RFs)){
  index <- index+1
  for(i in sequence(length(orderToGo)-1)){
    dataToUse <- which(orderToGo[i] == rownames(RFs))
    nextDataToUse <- which(orderToGo[i+1] == rownames(RFs))
    segments(orderedLevels[i], RFs[dataToUse, col], orderedLevels[i+1], RFs[nextDataToUse, col], col=cols[index])
  }
}
legend("topleft", legend=colnames(RFs), col=cols, lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
dev.off()



#  Figure 8 was made by hand


##  ----------------------------------------  ##
##             End Make Figures               ##
##  ----------------------------------------  ##




##  ----------------------------------------  ##
##            Supplement Figures              ##
##  ----------------------------------------  ##


#  Supplement Trees (14 datasets * 4 analyses)

setwd(FigDir)
pdf(file="SupplementTrees.pdf", width=8.5, height=11)
analyses <- unique(sapply(names(RAxML.TreeList), GetAnalysis, USE.NAMES=FALSE))
whichFig <- 0
for(whichAnal in analyses){
  red.trees <- RAxML.TreeList[grep(whichAnal, names(RAxML.TreeList))]
  if(whichAnal == "full")
    rr <- "All Sites"
  if(whichAnal == "nonasc")
    rr <- "Uncorrected"
  if(whichAnal == "lewis")
    rr <- "Conditional likelihood"
  if(whichAnal == "stam")
    rr <- "Reconstituted DNA"
  for(i in sequence(length(orderToGo))){
    whichFig <- whichFig+1
    nameToGrep <- paste0(whichAnal, "_out_", paste0(orderToGo[i], "noAmbigs"))
    if(whichAnal == "full")
      nameToGrep <- paste0(whichAnal, "_out_", paste0(orderToGo[i], "full"))
    tree <- RAxML.TreeList[[grep(nameToGrep, names(RAxML.TreeList))]]
    plot(tree, cex=0.5, edge.width=1, show.tip.label=TRUE, show.node.label=TRUE)
    add.scale.bar(cex = 0.7, font = 2, col = "black")
    bquote <- bquote(bold("Supplemental Figure S" * .(as.character(whichFig)) * ".  ") * .(as.character(rr)) * " RAxML tree for dataset " * .(as.character(orderToGo[i])) * ".")
    mtext(bquote, side=3, cex=0.75, adj=c(0))
  }
}
dev.off()

write.nexus(RAxML.TreeList, file="Trees.nex")



##  ----------------------------------------  ##
##         End Supplement Figures             ##
##  ----------------------------------------  ##




##  ----------------------------------------  ##
##      Double check RAxML Invocations        ##
##  ----------------------------------------  ##

# setwd(mainDir)
# systemCallSeeds <- CheckInvocations(getwd())
# CheckSeeds(systemCallSeeds)


##  ----------------------------------------  ##
##     End Double check RAxML Invocations     ##
##  ----------------------------------------  ##


