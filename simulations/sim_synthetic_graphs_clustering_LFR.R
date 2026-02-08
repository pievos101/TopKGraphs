library(TopKLists)
library(igraph)
library(fastcluster) 
library(cluster) 
library(netUtils)   # <-- For LFR benchmark graphs
library(aricode)
library(ggplot2)
library(reshape)
library(pheatmap)

# Load custom functions
source("/home/bpfeif/GitHub/TopKGraphs/simulations/sim.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/calc_BINARY.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs2.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bpfeif/GitHub/TopKGraphs/simulations/call_node2vec.R")

# Set seed for reproducibility
set.seed(123)

# ----------------------------
# Parameters for the simulation
# ----------------------------
n_iter <- 50    # Number of iterations
num_communities <- 3  # For evaluation purposes only (not used in LFR directly)

# Store results
RES_ari <- matrix(NaN, n_iter, 6)
colnames(RES_ari) <- c("TopKGraphs", "Jaccard", "Dice", "Laplacian", "PageRank", "Node2Vec")

RES_nmi <- matrix(NaN, n_iter, 6)
colnames(RES_nmi) <- c("TopKGraphs", "Jaccard", "Dice", "Laplacian", "PageRank", "Node2Vec")

# ----------------------------
# Simulation loop
# ----------------------------
for(xx in 1:n_iter){
  
  cat("\nIteration:", xx, "\n")
  
  # ----------------------------
  # Generate LFR benchmark graph in R
  # ----------------------------
  n_nodes <- 100
  avg_degree <- 5
  max_degree <- 10
  min_community <- 5
  max_community <- 50
  mu <- 0.05 # Mixing parameter
  tau1 <- 2 #2
  tau2 <- 1.1 #1.1   # tau2 must be > 1 for netUtils
  
  # Sample LFR graph
  g <- sample_lfr(n = n_nodes,
                  average_degree = avg_degree,
                  max_degree = max_degree,
                  min_community = min_community,
                  max_community = max_community,
                  mu = mu,
                  tau1 = tau1,
                  tau2 = tau2)
  

  print(table(V(g)$membership))

  # Ground truth communities
  V(g)$community <- V(g)$membership
  n_nodes <- vcount(g)
  
  # Remove isolated nodes (if any)
  g <- delete_vertices(g, V(g)[degree(g) == 0])
  n_nodes <- vcount(g)
  
  # ----------------------------
  # TopKGraphs
  # ----------------------------
  res <- topkgraphs(list(g), walk_depth=100, n_iter=50)
  
  topkgraphs_sim <- 1 - res$DIST
  
  # ----------------------------
  # Similarity matrices
  # ----------------------------
  jaccard_sim <- similarity(g, method = "jaccard", mode = "all")
  dice_sim <- similarity(g, method = "dice", mode = "all")
  
  # Laplacian embedding
  L <- laplacian_matrix(g, sparse = FALSE)
  eig <- eigen(L)
  emb <- eig$vectors
  
  # Node2Vec embedding
  node2vec_emb <- call_node2vec(g, walk_length=100, num_walks=50)

  node_order <- round(as.numeric(node2vec_emb[,1]))
  node2vec_emb <- as.matrix(node2vec_emb[,-1])
  ids <- match(1:n_nodes, node_order)
  node2vec_emb <- node2vec_emb[ids, ]
  
  # Personalized PageRank
  n <- vcount(g)
  ppr_mat <- sapply(seq_len(n), function(i) {
    pers <- rep(0, n)
    pers[i] <- 1
    page_rank(g, personalized = pers, damping=0.7)$vector
  })
  
  # ----------------------------
  # Hierarchical clustering
  # ----------------------------
  hc_topkgraphs <- hclust(as.dist(res$DIST), method="ward.D2")
  cl_topkgraphs <- cutree(hc_topkgraphs, num_communities)
  
  hc_jaccard <- hclust(as.dist(1 - jaccard_sim), method="ward.D2")
  cl_jaccard <- cutree(hc_jaccard, num_communities)
  
  hc_dice <- hclust(as.dist(1 - dice_sim), method="ward.D2")
  cl_dice <- cutree(hc_dice, num_communities)
  
  hc_laplacian <- hclust(dist(scale(emb)), method="ward.D2")
  cl_laplacian <- cutree(hc_laplacian, num_communities)
  
  hc_node2vec <- hclust(dist(scale(node2vec_emb)), method="ward.D2")
  cl_node2vec <- cutree(hc_node2vec, num_communities)
  
  hc_ppr <- hclust(as.dist(1 - ppr_mat), method="ward.D2")
  cl_ppr <- cutree(hc_ppr, num_communities)
  
  # ----------------------------
  # ARI and NMI evaluation
  # ----------------------------
  membership_gt <- V(g)$community
  
  ari_topkgraphs <- ARI(membership_gt, cl_topkgraphs)
  ari_jaccard   <- ARI(membership_gt, cl_jaccard)
  ari_dice      <- ARI(membership_gt, cl_dice)
  ari_laplacian <- ARI(membership_gt, cl_laplacian)
  ari_ppr       <- ARI(membership_gt, cl_ppr)
  ari_node2vec  <- ARI(membership_gt, cl_node2vec)
  
  nmi_topkgraphs <- NMI(membership_gt, cl_topkgraphs)
  nmi_jaccard   <- NMI(membership_gt, cl_jaccard)
  nmi_dice      <- NMI(membership_gt, cl_dice)
  nmi_laplacian <- NMI(membership_gt, cl_laplacian)
  nmi_ppr       <- NMI(membership_gt, cl_ppr)
  nmi_node2vec  <- NMI(membership_gt, cl_node2vec)
  
  # Store results
  RES_ari[xx,] <- c(ari_topkgraphs, ari_jaccard, ari_dice, ari_laplacian, ari_ppr, ari_node2vec)
  RES_nmi[xx,] <- c(nmi_topkgraphs, nmi_jaccard, nmi_dice, nmi_laplacian, nmi_ppr, nmi_node2vec)
  
  print(RES_ari)
}

# ----------------------------
# Plot results
# ----------------------------
RES <- RES_ari
RES_melt <- melt(RES)
colnames(RES_melt) <- c("ID","Method","Value")

p1 <- ggplot(RES_melt, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Clustering performance", y = "Adjusted R-Index (ARI)", x = "Method") +
  theme(legend.position = "none")

print(p1)

# ----------------------------
# Heatmap of TopKGraphs similarity
# ----------------------------
#pheatmap(topkgraphs_sim,
#         clustering_distance_rows = "euclidean",
#         clustering_distance_cols = "euclidean",
#         color = colorRampPalette(c("white", "blue"))(100),
#         main = "TopKGraphs Similarity Heatmap",
#         fontsize_row = 8,
#         fontsize_col = 8)

#cat("Simulation finished!\n")