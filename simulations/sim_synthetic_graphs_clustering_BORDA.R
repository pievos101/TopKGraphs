# Compare also to https://github.com/berenslab/graph-ne-paper
# Simulate synthetic graphs
# TopKGraphs
library(TopKLists)
library(igraph)
library(fastcluster) 
library(cluster) 
#library(node2vec)

source("/home/bastian/GitHub/TopKGraphs/simulations/sim.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_BINARY.R")
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bastian/GitHub/TopKGraphs/simulations/call_node2vec.R")

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
num_communities <- 3
nodes_per_community <- 10
intra_prob <- 0.8   # Probability of connection within communities
inter_prob <- 0.02  # Probability of connection between communities


n_iter = 50 

RES = matrix(NaN, n_iter, 4)
colnames(RES) = c("mean", "median", "geo.mean", "l2norm")

for(xx in 1:n_iter){

    #################################################
    # Generate the graph
    #g <- sample_islands(num_communities, 
    #                    nodes_per_community, 
    #                    5/10, 
    #                    2)
    
    # Plot the graph
    #plot(g, vertex.color = membership(cluster_label_prop(g)), 
    #    vertex.label = NA, vertex.size = 6, 
    #    main = "Graph with Known Community Structure")
    ######################################################

    # Sizes of the communities
    sizes <- c(10, 10, 10)
    n_nodes = sum(sizes)

    # Connection probability matrix (3x3)
    intra = 0.25 # 0.50 is baseline
    inter = 0.05 # 0.05 is baseline

    pref.matrix <- matrix(c(
      intra, inter, inter,
      inter, intra, inter,
      inter, inter, intra
    ), nrow = 3, byrow = TRUE)

    # Generate the SBM graph
    g <- sample_sbm(sum(sizes), pref.matrix, block.sizes = sizes)

    # Assign community membership
    V(g)$community <- rep(1:length(sizes), times = sizes)

    # Remove isolated nodes
    g <- delete_vertices(g, V(g)[degree(g) == 0])
    n_nodes = length(V(g))

    # Plot with communities
    #plot(g, vertex.color = V(g)$community, layout = layout_with_fr)

    res_mean = topkgraphs(list(g), walk_depth=20, n_iter=50, agg_method="mean")
    res_median = topkgraphs(list(g), walk_depth=20, n_iter=50, agg_method="median")
    res_geomean = topkgraphs(list(g), walk_depth=20, n_iter=50, agg_method="geo.mean")
    res_l2norm = topkgraphs(list(g), walk_depth=20, n_iter=50, agg_method="l2norm")
    

    hc_mean = hclust(as.dist(res_mean$DIST), method="ward.D2")
    cl_mean = cutree(hc_mean, num_communities)

    hc_median = hclust(as.dist(res_median$DIST), method="ward.D2")
    cl_median = cutree(hc_median, num_communities)

    hc_geomean = hclust(as.dist(res_geomean$DIST), method="ward.D2")
    cl_geomean = cutree(hc_geomean, num_communities)

    hc_l2norm = hclust(as.dist(res_l2norm$DIST), method="ward.D2")
    cl_l2norm = cutree(hc_l2norm, num_communities)

    
    library(aricode)
    #membership_gt <- rep(1:num_communities, each = nodes_per_community)
    membership_gt = V(g)$community  #rep(1:length(sizes), times = sizes)

    ari_mean = ARI(membership_gt, cl_mean)
    ari_median = ARI(membership_gt, cl_median)
    ari_geomean = ARI(membership_gt, cl_geomean)
    ari_l2norm = ARI(membership_gt, cl_l2norm)
    
    
    RES[xx,1] = ari_mean
    RES[xx,2] = ari_median
    RES[xx,3] = ari_geomean
    RES[xx,4] = ari_l2norm
  
    

print(RES)
}

library(ggplot2)
library(reshape)

RES_melt = melt(RES)
colnames(RES_melt) = c("ID","Method","Value")


p1 = ggplot(RES_melt, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Clustering performance",
    y = "Adjusted R-Index (ARI)",
    x = "Method"
  ) +
  theme(legend.position = "none")

print(p1)

stop("All good. finished!")

##################################################################
################### PLOTS ########################################








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

