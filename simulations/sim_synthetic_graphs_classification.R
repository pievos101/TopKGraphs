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
source("/home/bastian/GitHub/TopKGraphs/simulations/call_kNN.R")


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

RES = matrix(NaN, n_iter, 5)
colnames(RES) = c("TopKGraphs", "Jaccard", "Dice", "Laplacian", "Node2Vec")

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
    intra = 0.20 # 0.50 is baseline
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

    res = topkgraphs(list(g), walk_depth=5, n_iter=50)

    topkgraphs_sim = 1 - res$DIST

    jaccard_sim = similarity(g, method = "jaccard", mode = "all")

    dice_sim = similarity(g, method = "dice", mode = "all")

    # Laplace
    L    =  laplacian_matrix(g, sparse = FALSE)
    eig  = eigen(L)
    emb  = eig$vectors#[,1:10]

    #emb = embed_laplacian_matrix(g, 10)$X

    #print(emb)
    # Node2Vec
    #edges = get.edgelist(g)
    #emb = node2vecR(edges)
    #ids = as.numeric(rownames(emb))
    #ids = match(1:nrow(emb), ids)
    #emb = emb[ids, ]
    # ------------------

    # Call Python's Node2Vec
    node2vec_emb = call_node2vec(g, walk_length=5, num_walks=50)
    #print(dim(node2vec_emb))
    node_order = round(as.numeric(node2vec_emb[,1]))
    ids = match(1:length(V(g)), node_order)
    #print(node_order)
    node2vec_emb = as.matrix(node2vec_emb[ids,-1])

    #print(node2vec_emb)

    print(1)
    knn_topkgraphs = call_kNN_dist(res$DIST, V(g)$community, k = 5)
    print(2)
    knn_jaccard = call_kNN_dist(1-jaccard_sim,V(g)$community, k = 5)
    print(3)  
    knn_dice = call_kNN_dist(1-dice_sim,V(g)$community, k = 5)
    print(4)
    knn_laplacian = call_kNN_dist(as.matrix(dist(scale(emb))),
                                      V(g)$community, k = 5)
    print(5)
    knn_node2vec = call_kNN_dist(as.matrix(dist(scale(node2vec_emb))),
                                      V(g)$community, k = 5)
    print(6)
    # ----------------------------------------------- #
    #hc_node2vec = hclust(dist(emb), method="ward.D")
    #cl_node2vec = cutree(hc_node2vec, num_communities)
    #ids = as.numeric(names(cl_node2vec))
    #ids = match(1:n_nodes, node_order)
    #cl_node2vec = cl_node2vec[ids]
    #print(cl_node2vec)

    
    #library(aricode)
    #membership_gt <- rep(1:num_communities, each = nodes_per_community)
    membership_gt = rep(1:length(sizes), times = sizes)

    library(caret)
    all_levels <- sort(unique(c(knn_topkgraphs[[1]],knn_topkgraphs[[2]])))
    cm_topkgraphs = confusionMatrix(factor(knn_topkgraphs[[1]], levels=all_levels), 
                                    factor(knn_topkgraphs[[2]], levels=all_levels))

    all_levels <- sort(unique(c(knn_jaccard[[1]],knn_jaccard[[2]])))
    cm_jaccard = confusionMatrix(factor(knn_jaccard[[1]], levels=all_levels), 
                                    factor(knn_jaccard[[2]], levels=all_levels))

    all_levels <- sort(unique(c(knn_dice[[1]],knn_dice[[2]])))
    cm_dice = confusionMatrix(factor(knn_dice[[1]], levels=all_levels), 
                                    factor(knn_dice[[2]], levels=all_levels))

    all_levels <- sort(unique(c(knn_laplacian[[1]],knn_laplacian[[2]])))
    cm_laplacian = confusionMatrix(factor(knn_laplacian[[1]], levels=all_levels), 
                                    factor(knn_laplacian[[2]], levels=all_levels))

    all_levels <- sort(unique(c(knn_node2vec[[1]],knn_node2vec[[2]])))
    cm_node2vec = confusionMatrix(factor(knn_node2vec[[1]], levels=all_levels), 
                                    factor(knn_node2vec[[2]], levels=all_levels))
    
    RES[xx,1] = cm_topkgraphs$overall["Accuracy"]
    RES[xx,2] = cm_jaccard$overall["Accuracy"]
    RES[xx,3] = cm_dice$overall["Accuracy"]
    RES[xx,4] = cm_laplacian$overall["Accuracy"]
    RES[xx,5] = cm_node2vec$overall["Accuracy"]
    

print(RES)
}

stop("All good ....")

library(ggplot2)
library(reshape)

RES_melt = melt(RES)
colnames(RES_melt) = c("ID","Method","Value")


p1 = ggplot(RES_melt, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Classification performance",
    y = "Accuracy",
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

