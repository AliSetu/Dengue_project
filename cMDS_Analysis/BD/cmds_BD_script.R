# Load libraries
library(stats)
library(ape)
library(ade4)

# Load DNA sequences in the fasta format
dna <- read.dna(file = "DENV2_25_E_GENE_1485_594_aligned_BD_2006_2023.fasta", format = "fasta")
Dist <- dist.dna(dna, model = "TN93")

# Perform hierarchical clustering
hc <- hclust(Dist)
# Build phylogenetic tree (dendrogram) using hclust method
tree2 <- hclust(Dist)
# plot
plot(tree2, labels = NULL, hang = 0.1, check = TRUE,cex=0.6,  axes = TRUE, frame.plot = FALSE, ann = TRUE, main = "", sub = NULL, xlab = NULL, ylab = "Height")

# Cut the dendrogram to define clusters (adjust the height parameter as needed)
clusters <- cutree(hc, h = 0.015)  # Adjust height as needed

# Compute classical multidimensional scaling
mds <- cmdscale(Dist)

# Create a color palette for clusters
#length(unique(clusters))
cluster_colors <- rainbow(length(unique(clusters)))
j <- 0
# Plot the MDS results with points colored by clusters
plot(mds, type = "n", xlab = "Axis 1", ylab = "Axis 2", main = "Classical MDS of BD DNA sequences")
for (i in unique(clusters)) {
  points(mds[clusters == i, ], col = cluster_colors[i], pch = 19)
}

# Add legend
legend("bottomright", legend = unique(clusters), col = cluster_colors, pch = 19, title = "Cluster")

sink('cluster_bangladesh.txt')