##  ----------------------------------------  ##
##                                            ##
##           Phrynomics Analyses              ##
##           edited: 27 Oct 2014              ##
##                                            ##
##  ----------------------------------------  ##


library(phangorn)
library(phrynomics)
source("~/phrynomics-data/trunk/phrynomicsFunctions.R")


mainDir <- "/Users/Barb/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/newPhryno"
phrynoDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/newPhryno/vstam_vlewis_vfel_Results"
RepDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/RepeatRuns/"
FigDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/FULLfigs2"



##  ----------------------------------------  ##
##       Post-Analyses Trees and Data         ##
##  ----------------------------------------  ##

#  Load RAxML Trees and post-analyses scraping

setwd(phrynoDir)
files <- system("ls s*noAmbigs.phy", intern=T)
# analysis <- "RAxML"
# RAxML.trees <- system("ls RAxML_bipartitions.*", intern=T)  #trees with bootstraps
# RAxML.TreeList <- CreateTreeList(RAxML.trees, "RAxML")
# RAxML.TreeList <- lapply(RAxML.TreeList, ladderize)
# ML.results <- GetRAxMLStatsPostAnalysis(".")

# treeMatrix <- CreateTreeMatrix(RAxML.trees)
# ascgtr.treeMatrix <- AddTreeDist(treeMatrix[,-3], RAxML.TreeList)
# ascgtr.treeMatrix <- AddBLD(ascgtr.treeMatrix, RAxML.TreeList)
# ascgtr.BL.AllTrees <- list()
# for(i in sequence(dim(ascgtr.treeMatrix)[1])) {
  # tree1 <- assTrees(ascgtr.treeMatrix[i,1], RAxML.TreeList)[[1]]
  # tree2 <- assTrees(ascgtr.treeMatrix[i,2], RAxML.TreeList)[[1]]
  # ascgtr.BL.AllTrees[[i]] <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset=rownames(treeMatrix)[i])
  # names(ascgtr.BL.AllTrees)[[i]] <- rownames(ascgtr.treeMatrix)[i]
# }
# ascfull.treeMatrix <- AddTreeDist(treeMatrix[,c(3,1)], RAxML.TreeList)
# ascfull.treeMatrix <- AddBLD(ascfull.treeMatrix, RAxML.TreeList)
# ascfull.BL.AllTrees <- list()
# for(i in sequence(dim(ascfull.treeMatrix)[1])) {
  # tree1 <- assTrees(ascfull.treeMatrix[i,1], RAxML.TreeList)[[1]]
  # tree2 <- assTrees(ascfull.treeMatrix[i,2], RAxML.TreeList)[[1]]
  # ascfull.BL.AllTrees[[i]] <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset=rownames(treeMatrix)[i])
  # names(ascfull.BL.AllTrees)[[i]] <- rownames(ascfull.treeMatrix)[i]
# }


#load new vstam, vfel, vlewis
setwd(phrynoDir)
analysis <- "RAxML"
RAxML.trees <- system("ls RAxML_bipartitions.*", intern=T)  #trees with bootstraps
RAxML.TreeList <- CreateTreeList(RAxML.trees, "RAxML")
RAxML.TreeList <- lapply(RAxML.TreeList, ladderize)
ML.results <- GetRAxMLStatsPostAnalysis2(".")





#  Create global objects for tables and figs to be made

levels <-NULL
for(i in sequence(length(files))){
  levels <- as.numeric(c(levels, strsplit(files[i], "\\D+")[[1]][2])) #numerical datasets from file names
}
orderedLevels <- sort(levels)  #numerical datasets from file names
whichDatasets <- paste("s", levels, "", sep="")  #datasets from file names
AllOrder <- paste("s", seq(5, 65, 5), "", sep="") #datasets from sequence
orderToGo <- AllOrder[AllOrder %in% whichDatasets]  #datasets of sequence that exist
focalDatasets <- c("c5p3", "c25p3", "c55p3")


#setwd(mainDir)
#save(files, RAxML.trees, RAxML.TreeList, ML.results, ascgtr.treeMatrix, ascgtr.BL.AllTrees, ascfull.treeMatrix, ascfull.BL.AllTrees, MrBayes.trees, MrBayes.TreeList, MB.results, MB.ascgtr.treeMatrix, MB.ascgtr.BL.AllTrees, MB.ascfull.treeMatrix, MB.ascfull.BL.AllTrees, orderedLevels, orderToGo, focalDatasets, analyses, file="phrynoResults.Rdata")


# #for figures comparing ASC-GTR or FULL:ASC:
# comp="FULLASC"
# #comp <- "ASCGTR"
# if(comp == "ASCGTR"){
  # BL.AllTrees.RAxML <- ascgtr.BL.AllTrees 
  # BL.AllTrees.MrBayes <- MB.ascgtr.BL.AllTrees
  # treeMatrices <- list(ascgtr.treeMatrix, MB.ascgtr.treeMatrix)
# }
# if(comp == "FULLASC"){
  # BL.AllTrees.RAxML <- ascfull.BL.AllTrees
  # BL.AllTrees.MrBayes <- MB.ascfull.BL.AllTrees
  # treeMatrices <- list(ascfull.treeMatrix, MB.ascfull.treeMatrix)
# }


##  ----------------------------------------  ##
##      End Post-Analyses Trees and Data      ##
##  ----------------------------------------  ##






##  ----------------------------------------  ##
##                 Make Tables                ##
##  ----------------------------------------  ##

#  If you want to start from here, you can load up the workspace with all the data.
#load("phrynoResults.Rdata")

#  Table 1 was made by hand


#  Make Table 2. Summary ddRadSeq data. 

setwd(phrynoDir)
dataset <- sort(orderedLevels, decreasing=TRUE)
ASC.ML.results <- ML.results[which(ML.results[,2] == "ASC"),]
table4 <- matrix(nrow=length(dataset), ncol=6)
rownames(table4) <- paste("c", dataset, "p3", sep="")
colnames(table4) <- c("Matrix", "Loci", "VariableSites", "Missing(%)", "PhrynoOverlap", "non-PhrynoOverlap")
table4[,1] <- paste("s", dataset, sep="")
rows <- match(ASC.ML.results[,1], rownames(table4))
table4 <- table4[!is.na(rows),]  
table4[,c(2,3,4)] <- as.matrix(ASC.ML.results[order(rows[which(!is.na(rows))]), c(3,4,6)])
dataOverlap <- list()
for(i in sequence(length(files))){
  dataset <- read.table(files[[i]], row.names=1, colClasses="character", skip=1)
  dataOverlap[[i]] <- DataOverlap(dataset)[[2]]
  names(dataOverlap)[[i]] <- strsplit(strsplit(files[i], "/")[[1]][length(strsplit(files[i], "/")[[1]])], "no")[[1]][1]
}
for(m in sequence(dim(table4)[1])){
  whichDataset <- which(names(dataOverlap) == rownames(table4)[m])
  table4[m,6] <- round(mean(dataOverlap[[whichDataset]][-grep(pattern="PH", names(dataOverlap[[whichDataset]]))]), digits=2)
  table4[m,5] <- round(mean(dataOverlap[[whichDataset]][grep(pattern="PH", names(dataOverlap[[whichDataset]]))]), digits=2)
}
setwd(FigDir)
write.table(table4, file="table4.txt", quote=FALSE, row.names=FALSE)  


#  Table 5. Support for Clades

setwd(phrynoDir)
RAxTrees <- RAxML.TreeList[grep("c5p3|c25p3|c45p3|c65p3", names(RAxML.TreeList))]
MrBTrees <- MrBayes.TreeList[grep("c5p3|c25p3|c45p3|c65p3", names(MrBayes.TreeList))]
ascRAxTrees <- RAxTrees[grep("ASC", names(RAxTrees))]
gtrRAxTrees <- RAxTrees[grep("3_GTR", names(RAxTrees))]
fullRAxTrees <- RAxTrees[grep("full", names(RAxTrees))]
ascMrBTrees <- MrBTrees[grep("ASC", names(MrBTrees))]
gtrMrBTrees <- MrBTrees[grep("GTR_", names(MrBTrees))]
fullMrBTrees <- MrBTrees[grep("^c", names(MrBTrees))] 

Sceloporus <- c("SCMA1", "SCOC1", "SCGA1", "SCAN1")
Sceloporinae <- c(Sceloporus, "UROR1", "URBI1", "UTST1")
Tapaja <- c("PHHE5", "PHHE4", "PHHE2", "PHHE3", "PHHE1", "PHDI2", "PHDI1", "PHDO1", "PHDO2", "PHOR1", "PHOR2", "PHOR3", "PHOR4")
Doliosaurus <- c("PHMO3", "PHMO2", "PHGO3", "PHGO2", "PHGO1", "PHGO4", "PHPL1", "PHPL2", "PHPL3")
Brevicauda <- c("PHSP1", "PHSP4", "PHSP2", "PHSP3", "PHTA2", "PHTA3", "PHTA1", "PHTA4", "PHBR1", "PHBR3", "PHBR2", "PHBR4")
Anota <- c("PHMC1", "PHMC3", "PHMC2", "PHMC4", "PHSO3", "PHSO2", "PHSO1", "PHCE1", "PHCE3", "PHCE2", "PHCE4", "PHBL3", "PHBL2", "PHBL1", "PHBL4", "PHCE6", "PHCO1", "PHCE5")
Phrynosomatini <- c("PHAS4", "PHAS1", "PHAS3", "PHAS2", Anota, "PHCN1", "PHCN2", "PHCN3", "PHCN4", Brevicauda, Doliosaurus, Tapaja)
Callisaurini <- c("COTE1", "HOMA1", "CADR1", "CADR2")
Phrynosomatinae <- c(Callisaurini, Phrynosomatini)

clades <- list(Sceloporinae=Sceloporinae, Phrynosomatinae=Phrynosomatinae, Callisaurini=Callisaurini, Phrynosomatini=Phrynosomatini, Anota=Anota, Brevicauda=Brevicauda, Doliosaurus=Doliosaurus, Tapaja=Tapaja, Sceloporus=Sceloporus)

table5 <- matrix(nrow=length(clades), ncol=8)
rownames(table5) <- names(clades)
colnames(table5) <- paste0(rep(c("ML-", "BI-"), 4), c(rep("s65",2), rep("s45",2), rep("s25",2), rep("s5",2)))

for(i in sequence(dim(table5)[1])){
  taxa <- clades[[which(names(clades) == rownames(table5)[i])]]
  for(col in sequence(dim(table5)[2])){
    colnam <- colnames(table5)[col]
    level <- strsplit(colnam, "-")[[1]][2]
    level <- gsub("s", "c", level)
    threeVals <- rep("--", 3)
    if(length(grep("ML", colnam)) > 0){
      ascTr <- ascRAxTrees[[grep(level, names(ascRAxTrees))]]
      gtrTr <- gtrRAxTrees[[grep(level, names(gtrRAxTrees))]]
      fulTr <- fullRAxTrees[[grep(level, names(fullRAxTrees))]]
    }
    if(length(grep("BI", colnam)) > 0){
      ascTr <- ascMrBTrees[[grep(level, names(ascMrBTrees))]]
      gtrTr <- gtrMrBTrees[[grep(level, names(gtrMrBTrees))]]
      if(length(grep(level, names(fullMrBTrees))) > 0)  #some analyses not done
        fulTr <- fullMrBTrees[[grep(level, names(fullMrBTrees))]]
    }
    trees <- c(fulTr, gtrTr, ascTr)
    for(num in 1:3){
      tree <- trees[[num]]
      mrca <- getMRCA(tree, taxa)
      t2 <- nodeLeaves(tree, mrca)
      if(all(taxa %in% t2) && all(t2 %in% taxa)){
        EL <- GetEdgeList(tree)
        sup <- EL[which(EL[,2] == mrca), 5]
        threeVals[num] <- sup
      }
      if(length(grep("BI", colnam)) > 0){
        if(length(grep(level, names(fullMrBTrees))) == 0)
          threeVals[1] <- NA
      }
    }
    table5[i,col] <- paste(threeVals, collapse="/")
  }
}
setwd(FigDir)
write.table(table5, file="table5.txt", quote=FALSE, sep=" & ")
 

  

##  ----------------------------------------  ##
##               End Make Tables              ##
##  ----------------------------------------  ##






##  ----------------------------------------  ##
##                Make Figures                ##
##  ----------------------------------------  ##


#  Figure 1 was made by hand


#  Make Figure 2. Acquisition bias and tree length

setwd(FigDir)
pdf(file="Figure2.pdf", width=5, height=5)
layout(matrix(1:2, nrow=2, byrow=TRUE), respect=TRUE)
gtrcol <- "black"
asccol <- "gray"
fulcol <- "black"
gtrpoint <- "black"
ascpoint <- "gray"
fulpoint <- "white"
plot(rep(orderedLevels, 3), ML.results$TreeLength, type="n", ylab="Tree Length", xlab="Missing Data")
title(main="RAxML")
legend("topright", legend=c("Uncorrected", "Corrected", "All Sites"), col=c(gtrcol, asccol, fulcol), lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
for(i in sequence(length(orderToGo)-1)){
  dataToUse <- which(orderToGo[i] == ML.results$Level)
  nextDataToUse <- which(orderToGo[i+1] == ML.results$Level)
  segments(orderedLevels[i], ML.results$TreeLength[dataToUse[1]], orderedLevels[i+1], ML.results$TreeLength[nextDataToUse[1]], col=asccol)
  segments(orderedLevels[i], ML.results$TreeLength[dataToUse[2]], orderedLevels[i+1], ML.results$TreeLength[nextDataToUse[2]], col=gtrcol)  
  segments(orderedLevels[i], ML.results$TreeLength[dataToUse[3]], orderedLevels[i+1], ML.results$TreeLength[nextDataToUse[3]], col=fulcol)  
}
for(i in sequence(length(orderToGo))){
  dataToUse <- which(orderToGo[i] == ML.results$Level)
  points(orderedLevels[i], ML.results$TreeLength[dataToUse[1]], pch=21, bg=ascpoint)
  points(orderedLevels[i], ML.results$TreeLength[dataToUse[2]], pch=21, bg=gtrpoint)
  points(orderedLevels[i], ML.results$TreeLength[dataToUse[3]], pch=21, bg=fulpoint)
}
Ymin <- min(MB.results$TreeLength.lowCI, na.rm=TRUE)
Ymax <- max(MB.results$TreeLength.uppCI, na.rm=TRUE)
plot(rep(orderedLevels, 3), MB.results$TreeLength, type="n", ylab="Tree Length", xlab="Missing Data", ylim=c(Ymin, Ymax))
title(main="MrBayes")
legend("topright", legend=c("Uncorrected", "Corrected", "All Sites"), col=c(gtrcol, asccol, fulcol), lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
for(i in sequence(length(orderToGo)-1)){ 
  dataToUse <- which(orderToGo[i] == MB.results$Level)
  nextDataToUse <- which(orderToGo[i+1] == MB.results$Level)
  segments(orderedLevels[i], MB.results$TreeLength[dataToUse[1]], orderedLevels[i+1], MB.results$TreeLength[nextDataToUse[1]], col=asccol)
  segments(orderedLevels[i], MB.results$TreeLength[dataToUse[2]], orderedLevels[i+1], MB.results$TreeLength[nextDataToUse[2]], col=gtrcol)  
  segments(orderedLevels[i], MB.results$TreeLength[dataToUse[3]], orderedLevels[i+1], MB.results$TreeLength[nextDataToUse[3]], col=fulcol)  
}
for(i in sequence(length(orderToGo))){
  dataToUse <- which(orderToGo[i] == MB.results$Level)
  arrows(orderedLevels[i], MB.results$TreeLength.lowCI[dataToUse[1]], orderedLevels[i], MB.results$TreeLength.uppCI[dataToUse[1]], code=3, length=0.05, col=asccol, angle=90)
  arrows(orderedLevels[i], MB.results$TreeLength.lowCI[dataToUse[2]], orderedLevels[i], MB.results$TreeLength.uppCI[dataToUse[2]], code=3, length=0.05, col=gtrcol, angle=90)
  arrows(orderedLevels[i], MB.results$TreeLength.lowCI[dataToUse[3]], orderedLevels[i], MB.results$TreeLength.uppCI[dataToUse[3]], code=3, length=0.05, col=fulcol, angle=90)
  points(orderedLevels[i], MB.results$TreeLength[dataToUse[1]], pch=21, bg=ascpoint)
  points(orderedLevels[i], MB.results$TreeLength[dataToUse[2]], pch=21, bg=gtrpoint)
  points(orderedLevels[i], MB.results$TreeLength[dataToUse[3]], pch=21, bg=fulpoint)
}
dev.off()

#new Fig 2...AllSites, nonasc, vlewis, vstam, vfel
setwd(FigDir)
pdf(file="Figure2.pdf", width=5, height=8.5)
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



#  Make Figure 3. Scatterplot branch lengths

setwd(FigDir)
usr2dev <- function(x) 180/pi*atan(x*diff(par("usr")[1:2])/diff(par("usr")[3:4]))
getX <- function(y, linmod) (y-linmod$coefficients[1])/linmod$coefficients[2]  #(y -b)/m = x
pdf(file="Figure3.pdf", width=8.5, height=5)
op <- par(mar=par("mar")/1.7)
layout(matrix(1:8, nrow=2, byrow=TRUE), respect=TRUE)  
for(anal in 1:length(analyses)){
  whichAnalysis <- analyses[anal] 
  if(whichAnalysis == "RAxML")
    BL.AllTrees <- BL.AllTrees.RAxML
  if(whichAnalysis == "MrBayes")
    BL.AllTrees <- BL.AllTrees.MrBayes
  currmax <- 0
  for(i in sequence(length(focalDatasets))){
    dataToUse <- which(focalDatasets[i] == names(BL.AllTrees))
    if(length(dataToUse) > 0)
      currmax <- max(c(currmax, BL.AllTrees[[dataToUse]]$branchlength[which(BL.AllTrees[[dataToUse]]$present)], BL.AllTrees[[dataToUse]]$corr.BL[which(BL.AllTrees[[dataToUse]]$present)]))
  }
  for(i in sequence(length(focalDatasets))){
    dataToUse <- which(focalDatasets[i] == names(BL.AllTrees))
    if(length(dataToUse) == 0)
      plot(1:10, 1:10, type="n", axes=FALSE, frame.plot=FALSE, ylab="", xlab="")
    else {
      BLs <- BL.AllTrees[[dataToUse]]$branchlength[which(BL.AllTrees[[dataToUse]]$present)]
      corr.BLs <- BL.AllTrees[[dataToUse]]$corr.BL[which(BL.AllTrees[[dataToUse]]$present)]
      plot(BLs, corr.BLs, pch=21, bg="gray", ylim=c(0, currmax+currmax*.1), xlim=c(0, currmax+currmax*.1), xlab="ASC", ylab="non-ASC", type="n")
      linmod <- lm(corr.BLs ~ BLs)
      abline(linmod, lty=2)
      y <- 0.18
      points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
      points(BLs, corr.BLs, pch=21, bg="gray")
      text(x=getX(y, linmod), y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
      lines(c(-1,1), c(-1,1))
      title(main=paste("s", strsplit(focalDatasets[[i]], "\\D")[[1]][2], sep=""))
    }
  }
}
dev.off()


#  Figure 4 parts. Colored branch SNP phylogenies

setwd(FigDir)
figIndex <- 0
for(anal in 1:length(analyses)){
  whichAnalysis <- analyses[anal] 
  if(whichAnalysis == "RAxML"){
    TreeList <- RAxML.TreeList
    BL.AllTrees <- BL.AllTrees.RAxML
  }
  if(whichAnalysis == "MrBayes"){
    TreeList <- MrBayes.TreeList
    BL.AllTrees <- BL.AllTrees.MrBayes
  }
  for(i in sequence(length(orderToGo))){
    pdf(file=paste(whichAnalysis, ".", orderToGo[i], "trees.pdf", sep=""), width=8.5, height=11)
    figIndex <- figIndex+1
    dataToUse <- which(rownames(treeMatrices[[anal]]) == orderToGo[i])
    if(length(grep("NA", treeMatrices[[anal]][dataToUse,1])) > 0)
      plot(1:10, 1:10, type="n", axes=FALSE, frame.plot=FALSE, ylab="", xlab="")
    else {
      tree1 <- assTrees(treeMatrices[[anal]][dataToUse,1], TreeList)[[1]]
  print(is.rooted(tree1))
      tree2 <- assTrees(treeMatrices[[anal]][dataToUse,2], TreeList)[[1]]
      edgeColors <- BL.AllTrees[[dataToUse]]$edgeColor
      edgeColors[which(is.na(edgeColors))] <- rep("gray", length(which(is.na(edgeColors))))
    print(paste(names(BL.AllTrees)[i], mean(BL.AllTrees[[i]]$relativeBLdiff, na.rm=TRUE)))
      plot(tree1, edge.lty=BL.AllTrees[[dataToUse]]$edgelty, edge.color= edgeColors, cex=0.5, edge.width=2)
      legtxt <- c("Discordant", "< -10%", "-10% to 10%", "> 10%", "> 20%", "> 30%", "> 40%", "> 50%")
      legcolors <- c("gray", rgb(51,51,255, max=255), "gray", rgb(255,255,102, max=255), rgb(255,178,102, max=255), rgb(225,128,0, max=255), rgb(225,0,0, max=255), rgb(153,0,0, max=255))
      legend("bottomleft", legend=legtxt, col=legcolors, lty=c(2,rep(1,7)), lwd=2, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Branchlength Difference"))) 
      nodelabels(text=BL.AllTrees[[dataToUse]]$support[which(BL.AllTrees[[dataToUse]]$class == "internal")], node=BL.AllTrees[[dataToUse]]$desc[which(BL.AllTrees[[dataToUse]]$class == "internal")], cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, -0.1))
      nodelabels(text=BL.AllTrees[[dataToUse]]$corr.support[which(BL.AllTrees[[dataToUse]]$class == "internal")], node=BL.AllTrees[[dataToUse]]$desc[which(BL.AllTrees[[dataToUse]]$class == "internal")], cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, 1.1))
      if(whichAnalysis == "RAxML"){
        bquote1 <- bquote(bold("Supplemental Figure S" * .(as.character(figIndex)) * ".") * " Maximum likelihood phylogeny for the ascertainment-corrected dataset s" * .(strsplit(orderToGo[i], "\\D")[[1]][2]) ~ "(" * .(makeNumberwithcommas(ML.results[which(ML.results$Level == orderToGo[i]), 4][1])) ~ "variable sites," ~ "\n" )
        bquote2 <- bquote(.(unique(ML.results[which(ML.results$Level == orderToGo[i]), 6])) * "% missing data). Bootstrap values are shown on nodes (ascertainment-corrected above, no correction below). Branch")
        bquote3 <- bquote("color coding reflects the degree of relative length difference between ascertainment and non-ascertainment estimates.")
      }
      if(whichAnalysis == "MrBayes"){
        bquote1 <- bquote(bold("Supplemental Figure S" * .(as.character(figIndex)) * ".") * " Bayesian consensus phylogeny for the ascertainment-corrected dataset s" * .(strsplit(orderToGo[i], "\\D")[[1]][2]) ~ "(" * .(makeNumberwithcommas(MB.results[which(MB.results$Level == orderToGo[i]), 4][1])) ~ "variable sites," ~ "\n" )
        bquote2 <- bquote(.(unique(MB.results[which(MB.results$Level == orderToGo[i]), 6])) * "% missing data). Posterior probability values are shown on nodes (ascertainment-corrected above, no correction below). Branch")
        bquote3 <- bquote("color coding reflects the degree of relative length difference between ascertainment and non-ascertainment estimates.")
      }
      mtext(bquote1, side=3, cex=.75, adj=c(0), line=2)
      mtext(bquote2, side=3, cex=.75, adj=c(0), line=1)
      mtext(bquote3, side=3, cex=.75, adj=c(0), line=0)
  }
    dev.off()
  }
}


#  Make Figure 5. Mean branch length error (%)

setwd(FigDir)
pdf(file="Figure5.pdf", width=8.5, height=5)
for(anal in 1:length(analyses)){
  whichAnalysis <- analyses[anal] 
  if(whichAnalysis == "RAxML")
    BL.AllTrees <- BL.AllTrees.RAxML
  if(whichAnalysis == "MrBayes")
    BL.AllTrees <- BL.AllTrees.MrBayes
  taxon.BLdiff <- list()
  for(i in sequence(length(BL.AllTrees))){
    taxon.BLdiff[[i]] <- GetJustTipBLError(BL.AllTrees[[i]])
    names(taxon.BLdiff)[[i]] <- names(BL.AllTrees)[[i]]
  }
  mean.var.data <- NULL
  for(i in sequence(length(orderToGo))){
    whichData <- orderToGo[i]
    if(length(which(names(taxon.BLdiff) == whichData)) == 0)
      data4 <- rep(NA, 5)
    else {
      data1 <- taxon.BLdiff[[which(names(taxon.BLdiff) == whichData)]]
      data2 <- rep("OG", length(data1))
      data2[grep(pattern="PH", names(dataOverlap[[which(names(dataOverlap) == whichData)]]))] <- "PH"
      data3 <- data.frame(data1, data2)
      OGmean <- mean(data3[which(data3[,2] == "OG"),1])
      OGvar <- var(data3[which(data3[,2] == "OG"),1])
      PHmean <- mean(data3[which(data3[,2] == "PH"),1])
      PHvar <- var(data3[which(data3[,2] == "PH"),1])
      data4 <- c(whichData, OGmean, OGvar, PHmean, PHvar)
    }
    mean.var.data <- rbind(mean.var.data, data4)
  }
  colnames(mean.var.data) <- c("whichData", "OGmean", "OGvar", "PHmean", "PHvar")
  mean.var.data <- as.data.frame(mean.var.data, row.names= mean.var.data[,1], stringsAsFactors=FALSE)
  for(i in 2:5){
    mean.var.data[,i] <- as.numeric(mean.var.data[,i])
  }
  if(any(rownames(mean.var.data) == "c70p3"))
    mean.var.data <- mean.var.data[-which(rownames(mean.var.data) == "c70p3"),]
  if(whichAnalysis == "RAxML")
    RAxML.mean.var.data <- mean.var.data
  if(whichAnalysis == "MrBayes")
    MrBayes.mean.var.data <- mean.var.data
}
q <- sequence(dim(mean.var.data)[1])  #which(!is.na(mean.var.data[,1]))
par(mar=c(5,5,2,5))
plot(rep(q, 2), c(RAxML.mean.var.data[,2], MrBayes.mean.var.data[,2]), type="n", ylim=c(min(c(RAxML.mean.var.data[,4], RAxML.mean.var.data[,2], MrBayes.mean.var.data[,4], MrBayes.mean.var.data[,2]), na.rm=TRUE), max(c(RAxML.mean.var.data[,4], RAxML.mean.var.data[,2], MrBayes.mean.var.data[,4], MrBayes.mean.var.data[,2]), na.rm=TRUE)), ylab="mean", xlab="", axes=FALSE)
axis(side=2)
axis(side=1, at=q, labels=orderToGo[-which(orderToGo == "c70p3")])
points(q, RAxML.mean.var.data[,2], col="gray")
points(q, RAxML.mean.var.data[,4], col="black")
points(q, MrBayes.mean.var.data[,2], col="gray")
points(q, MrBayes.mean.var.data[,4], col="black")
for(i in sequence(dim(RAxML.mean.var.data)[1]-1)){
  segments(i, RAxML.mean.var.data[i,2], i+1, RAxML.mean.var.data[i+1,2], col="gray", lty=1)
  segments(i, RAxML.mean.var.data[i,4], i+1, RAxML.mean.var.data[i+1,4], col="black", lty=1)
}
for(i in sequence(dim(MrBayes.mean.var.data)[1]-1)){
  segments(i, MrBayes.mean.var.data[i,2], i+1, MrBayes.mean.var.data[i+1,2], col="gray", lty=2)
  segments(i, MrBayes.mean.var.data[i,4], i+1, MrBayes.mean.var.data[i+1,4], col="black", lty=2)
}
setwd(mainDir)
dev.off()


#  Make Figure 6 a-d. Scatterplot branch length support 
#  Figure 6 e-h was done using AWTY online to account for values < 50% and bipartition files not matching up correctly. 

setwd(FigDir)
pdf(file="Figure6a-d.pdf", width=8.5, height=5)
op <- par(mar=par("mar")/1.7)
layout(matrix(1:4, nrow=1, byrow=TRUE), respect=TRUE)
whichAnalysis <- "RAxML"
treeMatrix <- RAxML.treeMatrix
BL.AllTrees <- BL.AllTrees.RAxML
for(i in sequence(length(focalDatasets))){
  dataToUse <- which(focalDatasets[i] == names(BL.AllTrees))
  datarows <- which(BL.AllTrees[[dataToUse]]$present)[which(BL.AllTrees[[dataToUse]]$present) %in% which(BL.AllTrees[[dataToUse]]$support != 0)]  #don't want all the tip support 0s or non homologous clades
  support <- BL.AllTrees[[dataToUse]]$support[datarows]
  corr.support <- BL.AllTrees[[dataToUse]]$corr.support[datarows]
  xlims <- ylims <- c(1,100)
  plot(support, corr.support, ylab="GTR Tree", xlab="ASC Tree", pch=21, bg="gray", xlim=xlims, ylim=ylims)
  
  #text(support, corr.support, labels=BL.AllTrees[[dataToUse]][datarows,1])
  title(main=paste("s", strsplit(focalDatasets[[i]], "\\D")[[1]][2], sep=""))
}
setwd(mainDir)
dev.off()


#  Make Figure 7 a-b. Average RF distances among replicate runs.
#  Figure 7c was made using Mesquite. 

setwd(FigDir)
MLRFdistMatrix <- GetRFmatrix("RAxML")
BIRFdistMatrix <- GetRFmatrix("MrBayes")
RFdistMatrix <- list(ML=MLRFdistMatrix, BI=BIRFdistMatrix)
pdf(file="Figure7a-b.pdf", width=5, height=8.5)
layout(matrix(1:2, nrow=2, byrow=TRUE), respect=TRUE)
titles <- c("Ave. RF - 20 Independent Runs", "Ave. RF - 20 Posterior Trees")
for(dist in sequence(length(RFdistMatrix))){
  plot(as.numeric(RFdistMatrix[[dist]][,3]), as.numeric(RFdistMatrix[[dist]][,4]), xlab="dataset", ylab="Average RF", xaxt="n", type="n", ylim=c(0,130))
  title(main= titles[[dist]])
  cols <- c("lightblue", "blue")
  pchs <- c(15, 19)
  axis(side=1, at=1:14, labels=RFdistMatrix[[dist]][1:14, 2])
  ascsub <- RFdistMatrix[[dist]][which(RFdistMatrix[[dist]][,1] == "ASC"),]
  gtrsub <- RFdistMatrix[[dist]][which(RFdistMatrix[[dist]][,1] == "GTR"),]
  points(as.numeric(ascsub[,3]), as.numeric(ascsub[,4]), col=cols[1], pch=pchs[[dist]])
  points(as.numeric(gtrsub[,3]), as.numeric(gtrsub[,4]), col=cols[2], pch=pchs[[dist]])
  for(i in 1:13){
    segments(as.numeric(ascsub[i,3]), as.numeric(ascsub[i,4]), as.numeric(ascsub[i+1,3]), as.numeric(ascsub[i+1,4]), col=cols[1])
    segments(as.numeric(gtrsub[i,3]), as.numeric(gtrsub[i,4]), as.numeric(gtrsub[i+1,3]), as.numeric(gtrsub[i+1,4]), col=cols[2])
  }
  legtxt <- c("ASC", "GTR")
  legcolors <- c(cols[1], cols[2])
  legend("topleft", legend=legtxt, col=legcolors, lwd=2, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Model"))) 
}
dev.off()

#maybe add here how fig 7c was made? With scripts for exporting newick trees


#  Figure 8 was made by hand


##  ----------------------------------------  ##
##             End Make Figures               ##
##  ----------------------------------------  ##






##  ----------------------------------------  ##
##      Double check RAxML Invocations        ##
##  ----------------------------------------  ##

if(analysis == "RAxML"){
  systemCallSeeds <- CheckInvocations(getwd())
  CheckSeeds(systemCallSeeds)
}


##  ----------------------------------------  ##
##     End Double check RAxML Invocations     ##
##  ----------------------------------------  ##



## Old MrBayes stuff:

#  Load MrBayes Trees and post-analyses scraping  #not done YET...need to incorporate full too

analysis <- "MrBayes"
MrBayes.trees <- system(paste("ls *.con.tre", sep=""), intern=T)  
MrBayes.TreeList <- CreateTreeList(MrBayes.trees, "MrBayes")
MrBayes.TreeList <- lapply(MrBayes.TreeList, multi2di)  #remove later
MrBayes.TreeList <- lapply(MrBayes.TreeList, ladderize)  #remove later
MB.results <- GetMrBayesStatsPostAnalysis(".")  #will be warnings, because of missing analyses

treeMatrix <- CreateTreeMatrix(MrBayes.trees)
MB.ascgtr.treeMatrix <- AddTreeDist(treeMatrix[,-3], MrBayes.TreeList)
MB.ascgtr.treeMatrix <- AddBLD(MB.ascgtr.treeMatrix, MrBayes.TreeList)
MB.ascgtr.BL.AllTrees <- list()
for(i in sequence(dim(MB.ascgtr.treeMatrix)[1])) {
  tree1 <- assTrees(MB.ascgtr.treeMatrix[i,1], MrBayes.TreeList)[[1]]
  tree2 <- assTrees(MB.ascgtr.treeMatrix[i,2], MrBayes.TreeList)[[1]]
  MB.ascgtr.BL.AllTrees[[i]] <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset=rownames(treeMatrix)[i])
  names(MB.ascgtr.BL.AllTrees)[[i]] <- rownames(treeMatrix)[i]
}
MB.ascfull.treeMatrix <- AddTreeDist(treeMatrix[,c(3,1)], MrBayes.TreeList)
MB.ascfull.treeMatrix <- AddBLD(MB.ascfull.treeMatrix, MrBayes.TreeList)
MB.ascfull.BL.AllTrees <- list()
for(i in sequence(dim(MB.ascfull.treeMatrix)[1])) {
  if(length(grep("NA", MB.ascfull.treeMatrix[i,1])) == 0){
    tree1 <- assTrees(MB.ascfull.treeMatrix[i,1], MrBayes.TreeList)[[1]]
    tree2 <- assTrees(MB.ascfull.treeMatrix[i,2], MrBayes.TreeList)[[1]]
    MB.ascfull.BL.AllTrees[[i]] <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset=rownames(treeMatrix)[i])
    names(MB.ascfull.BL.AllTrees)[[i]] <- rownames(treeMatrix)[i]
  }
}
#MB.ascfull.BL.AllTrees <- MB.ascfull.BL.AllTrees[-c(which(names(MB.ascfull.BL.AllTrees) == "" ), which(is.na(names(MB.ascfull.BL.AllTrees))))]  #remove analyses that didn't finish (c5 and c10)






