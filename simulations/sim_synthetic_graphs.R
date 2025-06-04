# Simulate synthetic graphs
# TopKGraphs
library(TopKLists)
library(igraph)
library(fastcluster) 
library(cluster) 
library(node2vec)

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
num_communities <- 3
nodes_per_community <- 10
intra_prob <- 0.8   # Probability of connection within communities
inter_prob <- 0.02  # Probability of connection between communities


n_iter = 50 

RES = matrix(NaN, n_iter, 4)
colnames(RES) = c("TopKGraphs", "Jaccard", "Dice", "Laplacian")

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

    # Connection probability matrix (3x3)
    intra = 0.50
    inter = 0.25

    pref.matrix <- matrix(c(
      intra, inter, inter,
      inter, intra, inter,
      inter, inter, intra
    ), nrow = 3, byrow = TRUE)

    # Generate the SBM graph
    g <- sample_sbm(sum(sizes), pref.matrix, block.sizes = sizes)

    # Assign community membership
    V(g)$community <- rep(1:length(sizes), times = sizes)

    # Plot with communities
    #plot(g, vertex.color = V(g)$community, layout = layout_with_fr)

    res = topkgraphs(list(g), walk_depth=5, n_iter=20)

    topkgraphs_sim = 1 - res$DIST

    jaccard_sim = similarity(g, method = "jaccard", mode = "all")

    dice_sim = similarity(g, method = "dice", mode = "all")

    L    =  laplacian_matrix(g, sparse = FALSE)
    eig  = eigen(L)
    emb  = eig$vectors#[,1:10]

    # Node2Vec
    #edges = get.edgelist(g)
    #emb = node2vecR(edges)
    #ids = as.numeric(rownames(emb))
    #ids = match(1:nrow(emb), ids)
    #emb = emb[ids, ]

    hc_topkgraphs = hclust(as.dist(res$DIST), method="ward.D")
    cl_topkgraphs = cutree(hc_topkgraphs, num_communities)

    hc_jaccard = hclust(as.dist(1-jaccard_sim), method="ward.D")
    cl_jaccard = cutree(hc_jaccard, num_communities)

    hc_dice = hclust(as.dist(1-dice_sim), method="ward.D")
    cl_dice = cutree(hc_dice, num_communities)

    hc_laplacian = hclust(dist(scale(emb)), method="ward.D")
    cl_laplacian = cutree(hc_laplacian, num_communities)

    #hc_node2vec = hclust(dist(emb), method="ward.D")
    #cl_node2vec = cutree(hc_node2vec, num_communities)
    #ids = as.numeric(names(cl_node2vec))
    #ids = match(1:length(cl_node2vec), ids)
    #cl_node2vec = cl_node2vec[ids]
    #print(cl_node2vec)

    library(aricode)
    #membership_gt <- rep(1:num_communities, each = nodes_per_community)
    membership_gt = rep(1:length(sizes), times = sizes)

    ari_topkgraphs = ARI(membership_gt, cl_topkgraphs)
    ari_jaccard = ARI(membership_gt, cl_jaccard)
    ari_dice = ARI(membership_gt, cl_dice)
    ari_laplacian = ARI(membership_gt, cl_laplacian)
    
    #ari_node2vec = ARI(membership_gt, cl_node2vec)

    RES[xx,1] = ari_topkgraphs
    RES[xx,2] = ari_jaccard
    RES[xx,3] = ari_dice
    RES[xx,4] = ari_laplacian

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

