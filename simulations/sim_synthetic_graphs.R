# Simulate synthetic graphs
# TopKGraphs
library(TopKLists)
library(igraph)
library(fastcluster) 
library(cluster) 

source("/home/bastian/GitHub/TopKGraphs/simulations/sim.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_BINARY.R")
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs.R")

# Set seed for reproducibility
# set.seed(123)

# Parameters
#n_nodes <- 50      # Total number of nodes
#m_edges <- 2        # Number of edges to attach from a new node to existing nodes

# Generate Barabási–Albert graph
# g <- sample_pa(n = n_nodes, power = 1, m = m_edges, directed = FALSE)

# Plot the graph
#plot(g, vertex.size=5, vertex.label=NA, main="Barabási–Albert Graph")


# Parameters
num_communities <- 4
nodes_per_community <- 10
intra_prob <- 0.8   # Probability of connection within communities
inter_prob <- 0.02  # Probability of connection between communities

# Generate the graph
g <- sample_islands(num_communities, 
                    nodes_per_community, 
                    5/10, 
                    1)

# Plot the graph
plot(g, vertex.color = membership(cluster_label_prop(g)), 
     vertex.label = NA, vertex.size = 6, 
     main = "Graph with Known Community Structure")



res = topkgraphs(list(g), walk_depth=3, n_iter=20)
topkgraphs_sim = 1 - res$DIST

jaccard_sim <- similarity(g, method = "jaccard", mode = "all")

hc_topkgraphs = hclust(as.dist(res$DIST), method="ward.D")
cl_topkgraphs = cutree(hc_topkgraphs,4)

hc_jaccard = hclust(as.dist(1-jaccard_sim), method="ward.D")
cl_jaccard = cutree(hc_jaccard,4)

library(aricode)
membership_gt <- rep(1:num_communities, each = nodes_per_community)

ari_topkgraphs = ARI(membership_gt, cl_topkgraphs)
ari_jaccard = ARI(membership_gt, cl_jaccard)

print(ari_topkgraphs)
print(ari_jaccard)



################### PLOTS ###############################
# BASIC PLOT 
plot(hc)




# Plot heatmap
library(pheatmap)
pheatmap(topkgraphs_sim,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("white", "blue"))(100),
         main = "Similarity Heatmap",
         fontsize_row = 8,
         fontsize_col = 8)

