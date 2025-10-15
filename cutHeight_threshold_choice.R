
cutHeights <- seq(0.90, 0.999, by = 0.005)
module_summary <- data.frame(cutHeight = cutHeights, nModules = NA, medianSize = NA)

for (i in seq_along(cutHeights)) {
  ch <- cutHeights[i]
  
  dynamicMods <- cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    method = "hybrid",
    cutHeight = ch,
    deepSplit = TRUE,
    pamRespectsDendro = FALSE,
    minClusterSize = 5
  )
  
  # Assign colors and count module sizes
  colors <- labels2colors(dynamicMods)
  tab <- table(colors)
  
  module_summary$nModules[i] <- length(tab)
  module_summary$medianSize[i] <- median(tab)
}

# Visualize
par(mfrow = c(1, 2))
plot(module_summary$cutHeight, module_summary$nModules, type = "b",
     xlab = "cutHeight", ylab = "Number of modules", main = "Modules vs cutHeight")
plot(module_summary$cutHeight, module_summary$medianSize, type = "b",
     xlab = "cutHeight", ylab = "Median module size", main = "Median size vs cutHeight")
