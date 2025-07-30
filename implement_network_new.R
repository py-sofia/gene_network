################################### preprocessing ##############################


library(GEOquery)


expression_matrix <- read.csv("GSE128816.csv")
rownames(expression_matrix) <- expression_matrix[, 1]  # Set gene symbols as row names
expression_matrix <- expression_matrix[, -(1:2)]  # Remove first two columns
expression_matrix <- as.matrix(expression_matrix)  

expression_matrix <- log2(expression_matrix +1) # +1 to avoid log(0) (term: pseudocount)


summary(as.numeric(expression_matrix))

hist(as.numeric(expression_matrix), breaks = 100, main = "Value distribution")


# reduce complexity
# by excluding genes with low expression
avgexp <- rowMeans(expression_matrix, na.rm = TRUE)
expression_matrix <- expression_matrix[avgexp > quantile(avgexp, prob=0.25),] 
# keep genes with average expression above the 25th percentile





################################# GRN inference ################################



library(bc3net)
net.GSE128816 <- bc3net(expression_matrix, verbose=TRUE)


save(net.GSE128816, file="net.GSE128816.rda")


##
net.GSE128816 <- net.GSE4581
##

library(igraph)

net <- subgraph_from_edges(net.GSE128816,
                           eids = which(E(net.GSE128816)$weight>0.3))
net <- getgcc(net)

l <- layout_with_fr(net)
l <- norm_coords(l)

plot(net, layout=l, vertex.size=0.1,
     vertex.shape="none", vertex.label=NA, edge.width=0.6)

dcol = colorRampPalette(c("dodgerblue2", "white"))
points(l, col = densCols(l, colramp = dcol),
       pch = 20, cex=0.2,
       xaxt="n", yaxt="n", ylab="", xlab="")

getgcc = function(net) {
  mem <- which.max(clusters(net)$csize)
  genes.gcc <- V(net)$name[components(net)$membership == mem]
  net <- igraph::induced_subgraph(net, genes.gcc)
  
  return(net)
}
