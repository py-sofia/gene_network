

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


############################### weights histogram ############################## 

# Create data frame
df <- data.frame(
  weight = c(E(netTest)$weight, E(netControl)$weight),
  group = c(rep("netTest", ecount(netTest)),
            rep("netControl", ecount(netControl)))
)

# Plot
ggplot(df, aes(x = weight, fill = group)) +
  geom_histogram(position = "dodge", bins = 100, color = "white") +
  coord_cartesian(xlim = c(0, 0.25)) + 
  theme_minimal() +
  labs(title = "Edge Weight Distribution",
       x = "Edge Weight",
       y = "Count") +
  scale_fill_manual(values = c("skyblue", "tomato"))



################################ network graphs ################################ 

componentsInfoTest <- components(netTest)
giantTest <- induced_subgraph(netTest, 
                              which(componentsInfoTest$membership == which.max(componentsInfoTest$csize)))

comTest <- cluster_louvain(giantTest)
V(giantTest)$color <- brewer.pal(n = length(comTest), name = "Set3")[comTest$membership]

edgeColorFnTest <- col_numeric(palette = c("#fff", "#404040"), domain = NULL)
E(giantTest)$color <- edgeColorFnTest(E(giantTest)$weight)

pdf("netTest.pdf", width = 12, height = 12)
plot(giantTest,
     layout = layout_with_drl(giantTest),
     vertex.label = NA,
     vertex.size = 1,
     vertex.frame.color = NA,
     edge.arrow.mode = 0)
dev.off()



hist(E(netControl)$weight)

componentsInfoControl <- components(netControl)
giantControl <- induced_subgraph(netControl, 
                                 which(componentsInfoControl$membership == which.max(componentsInfoControl$csize)))

comControl <- cluster_louvain(giantControl)
V(giantControl)$color <- brewer.pal(n = length(comControl), name = "Set3")[comControl$membership]

edgeColorFnControl <- col_numeric(palette = c("#fff", "#404040"), domain = NULL)
E(giantControl)$color <- edgeColorFnControl(E(giantControl)$weight)

pdf("netzwerk.pdf", width = 12, height = 12)
plot(giantControl,
     layout = layout_with_drl(giantControl),
     vertex.label = NA,
     vertex.size = 1,
     vertex.frame.color = NA,
     edge.arrow.mode = 0)
dev.off()


################################ simple analysis ############################### 

# proportion of present edges from all possible edges in the network
edge_density(netTest, loops=F)
edge_density(netControl, loops=F)

# ratio of triangles to connected triples
transitivity(netTest, type="global")
transitivity(netControl, type="global")

# longest length of the shortest path between two nodes
diameter(netTest, directed=F) # ~2min to compute -> 2.04
diameter(netControl, directed=F) #               -> 1.73


################################### node degrees ###############################

# Degree vectors
degTest <- degree(netTest, mode = "all")
degControl <- degree(netControl, mode = "all")


# Combine into dataframe
df_deg <- data.frame(
  degree = c(degTest, degControl),
  group = c(rep("netTest", length(degTest)),
            rep("netControl", length(degControl)))
)

# Compute mean degrees
mean_df <- df_deg %>%
  group_by(group) %>%
  summarise(mean = mean(degree), .groups = "drop")

# Plot
ggplot(df_deg, aes(x = degree, fill = group)) +
  geom_histogram(
    binwidth = 1,
    color = "white",
    position = position_dodge(preserve = "single", width = 0.9)
  ) +
  coord_cartesian(xlim = c(0, 50)) +  # adjust as needed
  theme_minimal() +
  labs(
    title = "Node Degree Distribution (Side by Side)",
    x = "Degree",
    y = "Count"
  ) +
  scale_fill_manual(values = c("skyblue", "tomato")) +
  scale_color_manual(values = c("skyblue", "tomato"))


############################### degree distribution ############################

maxDeg <- 40

degControl.dist <- degree_distribution(netControl, cumulative=T, mode="all")
degTest.dist <- degree_distribution(netTest, cumulative=T, mode="all")

pdf("degree_distribution.pdf", width = 6, height = 4)
plot(x=0:maxDeg, y=1-degControl.dist[0:maxDeg+1], pch=19, cex=1.2, 
     col="skyblue", xlim=c(0, maxDeg), 
     xlab="Degree", ylab="Cumulative Frequency")
points(x=0:maxDeg, y=1-degTest.dist[0:maxDeg+1], pch=19, cex=1.2, col="tomato")
legend("bottomright", legend = c("Control", "Test"), col = c("skyblue", "tomato"), 
       pch = 19, inset=c(0, 0))
dev.off()



################################### centrality #################################






