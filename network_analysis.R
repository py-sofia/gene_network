# based on "DiffCoEx: a simple and sensitive method to find differentially coexpressed gene modules"






################################################################################
#                               data loading and preprocessing 
################################################################################



library(WGCNA) ### used for topological overlap calculation and clustering steps 
library(RColorBrewer) ### used to create nicer colour palettes 
library(preprocessCore) ### used by the quantile normalization function
library(igraph)
library(bc3net)
library(infotheo)
library(moduleColor)
library(flashClust)
library(ComplexHeatmap)


expressionMatrix <- read.csv("GSE128816.csv")

rownames(expressionMatrix) <- expressionMatrix[, 1]  # Set gene symbols as row names
expressionMatrix <- as.matrix(expressionMatrix[, 3:ncol(expressionMatrix)])  # Remove first two columns


normData <- normalize.quantiles(log2(expressionMatrix +1)) # +1 to avoid log(0) (term: pseudocount)

dimnames(normData) <- list(rownames(expressionMatrix), colnames(expressionMatrix))


datControl <- normData[, 1:10]  ### control group
datTreated <- normData[, 11:20] ### treatment group


################################################################################
#                             applying DiffCoEx
################################################################################




netControl <- bc3net(datControl, verbose=TRUE, estimator ="spearman", boot=100)
netTest <- bc3net(datTreated, verbose=TRUE, estimator="spearman", boot=100)


adjMatControl <- as_adjacency_matrix(netControl, attr="weight", sparse=F)
adjMatTreated <- as_adjacency_matrix(netTest, attr="weight", sparse=F)

genesTreated <- V(netTest)$name # get the names of all the treated genes
adjMatControl <- adjMatControl[genesTreated, genesTreated] # only keep genes present in treated

collectGarbage()


beta1=6 #user defined parameter for soft thresholding
# redo with pickSoftTreshold()

dissTOM <- TOMdist((abs(adjMatControl-adjMatTreated)/2)^(beta1/2))
rownames(dissTOM) <- rownames(adjMatControl)
colnames(dissTOM) <- rownames(adjMatControl)
collectGarbage()

#Hierarchical clustering is performed using the 
# Topological Overlap of the adjacency difference as input distance matrix
geneTree = flashClust(as.dist(dissTOM), method = "ward");

# Plot the resulting clustering tree (dendrogram)
png(file="hierarchicalTree.png",height=1000,width=1000)
plot(geneTreeC1C2, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels=F, hang = 0.04);
dev.off()

#We now extract modules from the hierarchical tree. This is done using cutreeDynamic. Please refer to WGCNA package documentation for details
dynamicModsHybrid = cutreeDynamic(dendro = geneTree, 
                                      distM = dissTOM,
                                      method="hybrid",
                                      cutHeight=0.999,
                                      deepSplit = T, 
                                      pamRespectsDendro = FALSE,
                                      minClusterSize = 5);

#Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
dynamicColorsHybrid = labels2colors(dynamicModsHybrid)

#the next step merges clusters which are close (see WGCNA package documentation)
mergedColor <- WGCNA::mergeCloseModules(t(cbind(datControl, datTreated)),
                                   dynamicColorsHybridC1C2,cutHeight=.2)$color
colorh1 <- mergedColor

#reassign better colors
colorh1[which(colorh1 =="midnightblue")]<-"red"
colorh1[which(colorh1 =="lightgreen")]<-"yellow"
colorh1[which(colorh1 =="cyan")]<-"orange"
colorh1[which(colorh1 =="lightcyan")]<-"green"

# Plot the dendrogram and colors underneath
png(file="module_assignment.png",width=1000,height=1000)
plotDendroAndColors(geneTree, 
                    colorh1, 
                    "Hybrid Tree Cut",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors cells")
dev.off()


# create anno (needed by extractModules())
geneSymbols <- rownames(expressionMatrix)
anno <- data.frame(gene_symbol = geneSymbols)
rownames(anno) <- geneSymbols


#We write each module to an individual file containing affymetrix probeset IDs
modulesMerged <- extractModules(colorh1$colors,
                                    t(datControl),
                                    anno,
                                    dir="modules",
                                    file_prefix=paste("Output","Specific_module",sep=''),
                                    write=T)
write.table(colorh1,
            file="module_assignment.txt",
            row.names=F,col.names=F,
            quote=F)

#We plot to a file the comparative heatmap showing correlation changes in the modules



################################################################################
#                               functions 
################################################################################

# uses module assignment list as input
# writes individual files witht he probeset ids for each module
extractModules <- 
  function(colorh1, datExpr, anno, write=F, file_prefix="", dir=NULL)
  {
    module <- list()
    if (!is.null(dir))
    {
      dir.create(dir)
      file_prefix=paste(dir, "/", file_prefix, sep="")
    }
    i <- 1
    for (c in unique(colorh1))
    {
      module[[i]] <- (anno[colnames(datExpr)[which(colorh1 == c)], 1])
      if (write)
      { 
        write.table(rownames(anno)[which(colorh1 == c)], 
                    file=paste(file_prefix, "_", c, ".txt", sep=""), 
                    quote=F, row.names=F, col.names=F)
      }
      i <- i+1
    }
    names(module) <- unique(colorh1)
    module
  }


# used by the plotting function to display close similar modules 
# based on their eigen values
getEigenGeneValues <- function(datRef, colorh1, datAll) 
{  
  eigenGenesCoef <- list() 
  i <- 0 
  for (c in unique(colorh1))
  {  
    i <- i+1 
    eigenGenesCoef[[i]] <- prcomp(
      scale(datRef[,which(colorh1 == c)]))$rotation[,1] 
  }
  names(eigenGenesCoef) <- unique(colorh1) 
  values <- NULL 
  for( c in unique(colorh1)) 
  {  
    v <- rbind(datAll)[,which(colorh1 == c)] %*% eigenGenesCoef[[c]] 
    values <- cbind(values,sign(mean(v))*v) 
  }  
  colnames(values) <- unique(colorh1) 
  values
}

# plotting function for comparative heatmap
plotC1C2Heatmap <- function(colorh1C1C2, AdjMat1C1, AdjMat1C2,datC1, datC2, 
                            ordering=NULL, file="DifferentialPlot.png")
{
  if (is.null(ordering))
  {
    h <- hclust(as.dist(1-abs(cor(getEigenGeneValues(
      datC1[,which(colorh1C1C2!="grey")], colorh1C1C2[which(colorh1C1C2!="grey")],
      rbind(datC1, datC2)[,which(colorh1C1C2!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering <- c(ordering, which(colorh1C1C2 ==c))
    }
  }
  mat_tmp <- (AdjMat1C1[ordering, ordering])
  mat_tmp[which(row(mat_tmp) > col(mat_tmp))] <- (AdjMat1C2[ordering,ordering]
                                                  [which(row(mat_tmp) > col(mat_tmp))])
  diag(mat_tmp) <- 0
  mat_tmp <- sign(mat_tmp)*abs(mat_tmp)^(1/2)
  png(file=file, height=1000, width=1000)
  image(mat_tmp, col=rev(brewer.pal(11,"RdYlBu")), axes=F, asp=1 ,breaks=
          seq(-1,1,length.out=12))
  dev.off()
  unique(colorh1C1C2[ordering])
}


# plots side by side the color bar of module assignments
# and the change in mean expression of the modules between the two conditions
plotExprChange<-function(datC1,datC2, colorhC1C2,ordering=NULL)
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,which(colorh1C1C2!="grey")],colorh1C1C2[which(colorh1C1C2!="grey")],rbind(datC1,datC2)[,which(colorh1C1C2!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1C1C2 ==c))
    }
  }
  mycolors<-colorh1C1C2[ordering]
  plot(x=0:length(which(mycolors!="grey")),y=rep(1,length(which(mycolors!="grey"))+1),col="white",axes=F,xlab="",ylab="",ylim=c(0,1))
  rr=c(244,239,225,215,209,193,181,166,151,130,110)
  gg=c(228,204,174,160,146,117,94,58,44,45,45)
  bb=c(176,140,109,105,102,91,84,74,70,68,66)
  MyColours<-NULL
  for ( i in 1:11)
  {
    MyColours=c(MyColours,rgb(rr[i],gg[i],bb[i],maxColorValue=255)  )
  }
  exprDiff<-NULL
  l<-0
  for (c in setdiff(unique(mycolors),"grey"))
  {
    meanC1<-mean(t(datC1)[colnames(datC1)[which(colorh1C1C2 == c)],])
    meanC2<-mean(t(datC2)[colnames(datC2)[which(colorh1C1C2 == c)],])
    exprDiff<-rbind(exprDiff,c(meanC1,meanC2))
    r<-l+length(which(mycolors==c))
    rect(l,0.85,r,1,col=c,border=F)
    rect(l,0,r,.4,col=MyColours[floor(meanC2*2)-10],border="white",lwd=2)
    rect(l,0.4,r,.8,col=MyColours[floor(meanC1*2)-10],border="white",lwd=2)
    l<-r
  }
  exprDiff
}


