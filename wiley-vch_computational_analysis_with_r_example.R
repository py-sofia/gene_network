
################################# IMPORTANT NOTES ##############################

# This code was adapted from the example provided in the Wiley-VCH book:
# Computational Network Analysis with R: Chapter 10

# some parts were modified to work with the latest versions of R packages

# the part "generating gene sets for functional analysis" was heavily changed 
# as the example code didn't work without converting the GO mappings from entrez 
# to gene symbols

# some snippets were added by me to improve speed
# I also added comments (with copilot)

################################ data preprocessing ############################

library(GEOquery)


eset <- getGEO('GSE4581', GSEMatrix = TRUE)[[1]] # retrieve data

data.probes <- exprs(eset) # extract expression values
data.probes <- log2(data.probes) # log2 transform the data
data.probes[which(is.na(data.probes))] <- 0 # replace NAs with 0


library(annotate)
library(hgu133plus2.db)

# either match probe set names to gene symbols and entrez gene identifiers
# entrez <- getEG(rownames(data.probes), data="hgu133plus2.db") # retrieve entrez gene identifiers
# entrez <- entrez[!is.na(entrez)] # remove NAs

# or assign gene symbols to probe set names
symbol <- getSYMBOL(rownames(data.probes), data="hgu133plus2.db") # retrieve gene symbols
data.probes <- data.probes[!is.na(symbol), ] # remove rows with NA symbols
symbol <- symbol[!is.na(symbol)] # remove NAs from symbols


# gene level summarization
names(symbol) <- rownames(data.probes)
set <- split(names(symbol), symbol)

expmat <- matrix(NA, nrow=length(set), ncol=ncol(data.probes)) # create an empty matrix
colnames(expmat) <- colnames(data.probes) # assign column names
rownames(expmat) <- names(set) # assign row names

for (i in 1:length(set)) {
  cat(i, "\n")
  probes <- set[[i]] # get the probe set for the i-th gene
  if (length(probes) == 1) { # if there is only one probe for the gene
    expmat[i, ] <- data.probes[probes, ] # assign the expression value of that probe
  } else {
    expmat[i, ] <- apply(data.probes[probes,], 2, median) # compute the median expression value across probes
  }
}

save(expmat, file="expmat.rda")


# DOESN'T WORK
# reduce complexity
# by excluding genes with low expression
# avgexp <- rowMeans(expmat, na.rm = TRUE)
# expmat <- expmat[avgexp > quantile(avgexp, prob=0.25),] 
        # keep genes with average expression above the 25th percentile

# Keep top 1000 most highly expressed genes
avgexp <- rowMeans(expmat, na.rm = TRUE)
top_genes <- order(avgexp, decreasing = TRUE)[1:1000]
expmat <- expmat[top_genes,]


################################# GRN inference ################################

library(bc3net)
net.GSE4581 <- bc3net(expmat, verbose=TRUE)

save(net.GSE4581, file="net.GSE4581.rda")

library(igraph)

net <- subgraph_from_edges(net.GSE4581,
                      eids = which(E(net.GSE4581)$weight>0.3))
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



############## generating gene sets for functional analysis ####################

library("org.Hs.eg.db")
library("GO.db")

# 1. Retrieve GO annotations (Entrez IDs by default)
GO.frame <- toTable(org.Hs.egGO2ALLEGS) # gene_id (Entrez) to GO mapping
GO.frame <- GO.frame[!duplicated(GO.frame), ] # remove duplicates

# 2. Map Entrez IDs to Gene Symbols
entrez2symbol <- toTable(org.Hs.egSYMBOL)  # Entrez → Symbol mapping
colnames(entrez2symbol) <- c("gene_id", "symbol")

GO.frame <- merge(GO.frame, entrez2symbol, by = "gene_id", all.x = TRUE)
GO.frame <- GO.frame[!is.na(GO.frame$symbol), ] # keep only rows with symbols

# 3. Build GO term list using gene symbols
GO <- split(GO.frame$symbol, GO.frame$go_id)

# 4. Add descriptive names for GO terms
goterm <- mget(names(GO), GOTERM, ifnotfound = NA)
term <- sapply(goterm, function(x) if (isS4(x)) x@Term else "")
names(GO) <- paste(names(GO), ":", term)

# 5. Check overlap between network genes and GO gene sets
network_genes <- V(net.GSE4581)$name
overlap_count <- length(intersect(network_genes, unique(unlist(GO))))
cat("Overlap between network genes and GO sets:", overlap_count, "\n")



########################### pathway and other gene set collections #############
# skipped


######################## functional enrichment analysis ########################

res <- gpea(net.GSE4581, GO, cmax=1000, cmin=2) # perform GSEA on the network

save(res, file="res.GSE4581.rda")


#################### visualization of single network component #################

net <- subgraph_from_edges(net.GSE4581,
                      eids = which(E(net.GSE4581)$weight>0.15))


genes <- GO[["GO:0048583 : regulation of response to stimulus"]]

present <- V(net)$name%in%genes

module <- getgcc(induced_subgraph(net, V(net)$name[present]))

plot(module, layout=layout.fruchterman.reingold, vertex.size=2,
     vertex.label.cex=0.35,
     vertex.label.color="black",
     vertex.color="dodgerblue",
     vertex.frame.color="dodgerblue2",
     vertex.label.dist=0.2)








