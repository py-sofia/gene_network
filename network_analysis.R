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
library(scales)
library(ggplot2)
library(dplyr)


expressionMatrix <- read.csv("GSE128816.csv")

rownames(expressionMatrix) <- expressionMatrix[, 1]  # Set gene symbols as row names
expressionMatrix <- as.matrix(expressionMatrix[, 3:ncol(expressionMatrix)])  # Remove first two columns


normData <- normalize.quantiles(log2(expressionMatrix +1)) # +1 to avoid log(0) (term: pseudocount)

dimnames(normData) <- list(rownames(expressionMatrix), colnames(expressionMatrix))


datControl <- normData[, 1:10]  ### control group
datTreated <- normData[, 11:20] ### treatment group






################################################################################
#                             applying DiffCoEx (part 1)
################################################################################



netControl <- bc3net(datControl, verbose=TRUE, estimator="spearman", boot=100)
netTest <- bc3net(datTreated, verbose=TRUE, estimator="spearman", boot=100)
# ca. 24h to run 

save(netControl, file = "netControl.RData")
save(netTest, file = "netTest.RData")


adjMatControl <- as_adjacency_matrix(netControl, attr="weight", sparse=F)
adjMatTreated <- as_adjacency_matrix(netTest, attr="weight", sparse=F)

genesTreated <- V(netTest)$name # get the names of all the treated genes
adjMatControl <- adjMatControl[genesTreated, genesTreated] # only keep genes present in treated

collectGarbage()

################################################################################
#                             pickSoftThreshold()
################################################################################


AdjDiff <- abs(adjMatControl-adjMatTreated)/2
diag(AdjDiff) <- 1

powers <- c(seq(1,10,1), seq(12,20,2))

sft = pickSoftThreshold.fromSimilarity(
  similarity = AdjDiff,
  powerVector = powers,
  RsquaredCut = 0.85,
  moreNetworkConcepts = TRUE,
  verbose = 5);

beta1 <- sft$powerEstimate # optimal beta1
beta1 <- 3

head(sft$fitIndices)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.85, col="red")



################################################################################
#                             applying DiffCoEx (part 2)
################################################################################


dissTOM <- TOMdist((AdjDiff)^(beta1/2))
rownames(dissTOM) <- rownames(adjMatControl)
colnames(dissTOM) <- rownames(adjMatControl)
collectGarbage()

#Hierarchical clustering is performed using the 
# Topological Overlap of the adjacency difference as input distance matrix
geneTree = flashClust(as.dist(dissTOM), method = "ward");

# Plot the resulting clustering tree (dendrogram)
png(file="hierarchicalTree.png",height=1000,width=3000)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels=F, hang = 0.04);
dev.off()

#We now extract modules from the hierarchical tree. This is done using cutreeDynamic. Please refer to WGCNA package documentation for details
dynamicModsHybrid = cutreeDynamic(dendro = geneTree, 
                                      distM = dissTOM,
                                      method="hybrid",
                                      cutHeight=0.96,
                                      deepSplit = T, 
                                      pamRespectsDendro = FALSE,
                                      minClusterSize = 5);

#Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
dynamicColorsHybrid = labels2colors(dynamicModsHybrid)

#the next step merges clusters which are close (see WGCNA package documentation)
mergedColor <- WGCNA::mergeCloseModules(t(cbind(datControl, datTreated)),
                                   dynamicColorsHybrid,cutHeight=.2)$color
colorh1 <- mergedColor

#reassign better colors
colorh1[which(colorh1 =="midnightblue")]<-"red"
colorh1[which(colorh1 =="lightgreen")]<-"yellow"
colorh1[which(colorh1 =="cyan")]<-"orange"
colorh1[which(colorh1 =="lightcyan")]<-"green"

# Plot the dendrogram and colors underneath
png(file="module_assignment.png",width=4000,height=1000)
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
modulesMerged <- extractModules(colorh1,
                                    t(datControl),
                                    anno,
                                    dir="modules",
                                    file_prefix=paste("Output","Specific_module",sep=''),
                                    write=T)

write.table(colorh1,
            file="module_assignment.txt",
            row.names=F,col.names=F,
            quote=F)



################################################################################
#                             comparative heatmap
################################################################################

library(magick)
library(circlize)

col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

pdf("AdjDiff_heatmap.pdf", width = 8, height = 8)
Heatmap(AdjDiff[1:2000, 1:2000], col=col_fun,
        show_row_dend=F,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=F,
        cluster_rows=T,
        cluster_columns=T,
        use_raster=T,
        raster_quality=1)
dev.off()

