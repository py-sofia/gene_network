

library(clusterProfiler)
library(org.Hs.eg.db)  # assuming human genes


modulesMergedEntrez <- lapply(modulesMerged, function(gene_set) {
  result <- bitr(gene_set,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb = org.Hs.eg.db)
  return(result$ENTREZID)
})


enrichment_results <- lapply(modulesMergedEntrez, function(gene_set) {
  enrichGO(gene=modulesMergedEntrez$gene_set,
           OrgDb=org.Hs.eg.db,
           ont="BP",
           pAdjustMethod="BH",
           readable=T,
           pvalueCutoff=0.05)
})
