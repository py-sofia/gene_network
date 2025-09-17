# based on "DiffCoEx: a simple and sensitive method to find differentially coexpressed gene modules"


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



################################################################################
#                               data loading and preprocessing 
################################################################################



library(WGCNA) ### used for topological overlap calculation and clustering steps 
library(RColorBrewer) ### used to create nicer colour palettes 
library(preprocessCore) ### used by the quantile normalization function


expression_matrix <- read.csv("GSE128816.csv")
rownames(expression_matrix) <- expression_matrix[, 1]  # Set gene symbols as row names


expression_matrix <- expression_matrix[, -(1:2)]  # Remove first two columns
expression_matrix <- as.matrix(expression_matrix)  


norm_data <- normalize.quantiles(log2(expression_matrix +1)) # +1 to avoid log(0) (term: pseudocount)
dimnames(norm_data) <- list(rownames(expression_matrix), colnames(expression_matrix))


datC1 <- norm_data[, 1:10]  ### control group
datC2 <- norm_data[, 11:20] ### treatment group




################################################################################
#                             applying DiffCoEx
################################################################################

beta1=6 #user defined parameter for soft thresholding

# redo with pickSoftTreshold()

# import from BC3NET
library(igraph)
library(bc3net)
library(infotheo)
library(moduleColor)

net_control <- bc3net(datC1, verbose=TRUE, estimator="emp")
net_test <- bc3net(datC2, verbose=TRUE, estimator="emp")


#load("~/Schule/Matura/MA/Modellierung/gene_network/net_control_group.rda")
AdjMatC1 <- as_adjacency_matrix(net_control, attr="weight", sparse=F)

#load("~/Schule/Matura/MA/Modellierung/gene_network/net_test_group.rda")
AdjMatC2 <- as_adjacency_matrix(net_test, attr="weight", sparse=F)


# Extract vertex names
genesC1 <- V(net_control)$name
genesC2 <- V(net_test)$name

# Reorder AdjMatC2 to match AdjMatC1
AdjMatC2 <- AdjMatC2[genesC1, genesC1]



diag(AdjMatC1)<-0
diag(AdjMatC2)<-0
collectGarbage()

dissTOMC1C2=TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(beta1/2))
collectGarbage()

#Hierarchical clustering is performed using the 
# Topological Overlap of the adjacency difference as input distance matrix
library(flashClust)
geneTreeC1C2 = flashClust(as.dist(dissTOMC1C2), method = "average");

# Plot the resulting clustering tree (dendrogram)
png(file="hierarchicalTree.png",height=1000,width=1000)
plot(geneTreeC1C2, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
dev.off()

#We now extract modules from the hierarchical tree. This is done using cutreeDynamic. Please refer to WGCNA package documentation for details
################################################# parameters are critical (esp. cutHeight, minClusterSize)
dynamicModsHybridC1C2 = cutreeDynamic(dendro = geneTreeC1C2, distM = dissTOMC1C2,method="hybrid",cutHeight=0.999,deepSplit = T, pamRespectsDendro = FALSE,minClusterSize = 5);

#Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
dynamicColorsHybridC1C2 = labels2colors(dynamicModsHybridC1C2)

#the next step merges clusters which are close (see WGCNA package documentation)
################################################# parameters are critical (esp. cutHeight)
mergedColorC1C2<-mergeCloseModules(t(norm_data),dynamicColorsHybridC1C2,cutHeight=.1)$color
colorh1C1C2<-mergedColorC1C2

#reassign better colors
colorh1C1C2[which(colorh1C1C2 =="midnightblue")]<-"red"
colorh1C1C2[which(colorh1C1C2 =="lightgreen")]<-"yellow"
colorh1C1C2[which(colorh1C1C2 =="cyan")]<-"orange"
colorh1C1C2[which(colorh1C1C2 =="lightcyan")]<-"green"
# Plot the dendrogram and colors underneath
png(file="module_assignment.png",width=1000,height=1000)
plotDendroAndColors(geneTreeC1C2, colorh1C1C2, "Hybrid Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors cells")
dev.off()


# create anno (needed by extractModules())
gene_symbols <- rownames(expression_matrix)
anno <- data.frame(gene_symbol = gene_symbols)
rownames(anno) <- gene_symbols


#We write each module to an individual file containing affymetrix probeset IDs
modulesC1C2Merged<-extractModules(colorh1C1C2,t(datC1),anno,dir="modules",file_prefix=paste("Output","Specific_module",sep=''),write=T)
write.table(colorh1C1C2,file="module_assignment.txt",row.names=F,col.names=F,quote=F)

#We plot to a file the comparative heatmap showing correlation changes in the modules
#The code for the function plotC1C2Heatmap and others can be found below under the Supporting Functions section

plotC1C2Heatmap(colorh1C1C2,AdjMatC1,AdjMatC2, t(datC1), t(datC2))
png(file="exprChange.png",height=500,width=500)
plotExprChange(t(datC1),t(datC2),colorh1C1C2)
dev.off()




