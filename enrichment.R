

library(clusterProfiler)
library(org.Hs.eg.db)  # assuming human genes
library(enrichplot)



modulesMergedEntrez <- lapply(modulesMerged, function(gene_set) {
  result <- bitr(gene_set,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb = org.Hs.eg.db)
  return(result$ENTREZID)
})



universe <- bitr(geneSymbols,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb = org.Hs.eg.db)

enrichmentResults_ontALL <- lapply(modulesMergedEntrez, function(gene_set) {
  enrichGO(
    gene = gene_set,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    readable = TRUE,
    pvalueCutoff = 0.05,
    qvalueCutoff = 1,
    universe = universe
  )
})



lapply(names(enrichmentResults_ontALL), function(module_name) {
  enrich_result <- enrichmentResults_ontALL[[module_name]]
  
  if (nrow(as.data.frame(enrich_result)) > 0) {  # only plot if non-empty
    print(barplot(enrich_result, showCategory = 20, title = module_name))
    Sys.sleep(5)
  } else {
    message(paste("No enrichment results for module:", module_name))
  }
})



compareClusterResult_ontALL <- compareCluster(
  geneClusters=modulesMergedEntrez,
  fun="enrichGO",
  OrgDb=org.Hs.eg.db,
  keyType="ENTREZID",
  ont="ALL",
  pAdjustMethod="BH", # p-values adjusted by cluster
  readable=T,
  pvalueCutoff=0.05,
  qvalueCutoff=1, # do not filter by q-value
  universe=universe) 

barplot(compareClusterResult_ontALL, showCategory=20)
