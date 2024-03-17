# Load libraries
library(stats)
library(ape)
library(ade4)

# Load DNA sequences in the fasta format
dna <- read.dna(file = "global_Cosmopolitan_C_aligned.fasta", format = "fasta")
Dist <- dist.dna(dna, model = "TN93")

# Perform hierarchical clustering
hc <- hclust(Dist)
# Build phylogenetic tree (dendrogram) using hclust method
tree2 <- hclust(Dist)
# plot
plot(tree2, labels = NULL, hang = 0.1, check = TRUE,cex=0.6,  axes = TRUE, frame.plot = FALSE, ann = TRUE, main = "", sub = NULL, xlab = NULL, ylab = "Height")

# Cut the dendrogram to define clusters (adjust the height parameter as needed)
clusters <- cutree(hc, h = 0.025)  # Adjust height as needed

# Compute classical multidimensional scaling
mds <- cmdscale(Dist)

# Create a color palette for clusters
#length(unique(clusters))
cluster_colors <- rainbow(length(unique(clusters)))
j <- 0
# Plot the MDS results with points colored by clusters
plot(mds, type = "n", xlab = "Axis 1", ylab = "Axis 2", main = "Classical MDS of DNA sequences")
for (i in unique(clusters)) {
  points(mds[clusters == i, ], col = cluster_colors[i], pch = 19)
}

# Add legend
#legend("bottomright", legend = unique(clusters), col = cluster_colors, pch = 19, title = "Cluster")
#legend("bottomright", 
       #legend = unique(clusters), 
       #col = cluster_colors, 
       #pch = 19, 
       #title = "Cluster",
       #inset=c(-0.1,0))
# Increase margin space for the legend
par(xpd=TRUE, mar = c(8, 8, 8, 8))

# Add legend outside the plot
#legend("bottomright", 
       #legend = unique(clusters), 
       #col = cluster_colors, 
       #pch = 19, 
       #title = "Cluster",
       #inset=c(0,0))
plot(x_values, y_values)

# Use locator to click on the plot and obtain coordinates
legend_coords <- locator(1)

# Add legend using the coordinates obtained from locator
legend(legend_coords$x, legend_coords$y,  # Use coordinates obtained from locator
       legend = unique(clusters), 
       col = cluster_colors, 
       pch = 19, 
       title = "Clusters", 
       cex = .6,
       bg = "white",  # Optional: set background color of the legend
       box.lwd = 1)   # Optional: set border width of the legend box
