# Compare also to https://github.com/berenslab/graph-ne-paper
# Simulate synthetic graphs
# TopKGraphs
library(TopKLists)
library(igraph)
library(fastcluster) 
library(cluster) 
#library(node2vec)

source("/home/bpfeif/GitHub/TopKGraphs/simulations/sim.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/calc_BINARY.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bpfeif/GitHub/TopKGraphs/simulations/call_node2vec.R")

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

RES = matrix(NaN, n_iter, 6)
colnames(RES) = c("TopKGraphs", "Jaccard", "Dice", "Laplacian", "PageRank",
                "Node2Vec")

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
    sizes = c(20,20,20) #sample(c(5,10,20),3,replace=TRUE) #c(10, 20, 5)
    n_nodes = sum(sizes)

    # Connection probability matrix (3x3)
    intra = 0.99 # 0.50 is baseline
    inter = 0.70 # 0.05 is baseline

    pref.matrix <- matrix(c(
      intra, inter, inter,
      inter, intra, inter,
      inter, inter, intra
    ), nrow = 3, byrow = TRUE)

    # Generate the SBM graph
    g <- sample_sbm(sum(sizes), pref.matrix, block.sizes = sizes)

    # Assign community membership
    V(g)$community <- rep(1:length(sizes), times = sizes)

    ### optional - perturbation
    # Remove fraction p of edges
    #p <- 0.4
    #edges_to_remove <- E(g)[runif(ecount(g)) < p]
    #g <- delete_edges(g, edges_to_remove)
    
    # Rewire 10% of edges
    #frac <- 0.3
    #num_rewire <- ceiling(frac * ecount(g))
    #g <- rewire(g, with = keeping_degseq(niter = num_rewire))
    ####

    # Remove isolated nodes
    g = delete_vertices(g, V(g)[degree(g) == 0])
    n_nodes = length(V(g))

    # Plot with communities
    #plot(g, vertex.color = V(g)$community, layout = layout_with_fr)

    res = topkgraphs(list(g), walk_depth=20, n_iter=50, 
                          do.BORDA = TRUE,
                          do.TopKSignal=FALSE,
                          do.RRA = FALSE)
    #print(res)

    topkgraphs_sim = 1 - res$DIST

    jaccard_sim = similarity(g, method = "jaccard", mode = "all")

    dice_sim = similarity(g, method = "dice", mode = "all")

    # Laplace
    L    = laplacian_matrix(g, sparse = FALSE)
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
    node2vec_emb = call_node2vec(g, walk_length=20, num_walks=50)
    #print(dim(node2vec_emb))
    node_order = round(as.numeric(node2vec_emb[,1]))
    #print(node_order)
    node2vec_emb = as.matrix(node2vec_emb[,-1])

    #print(node2vec_emb)

    # Page_rank
    n <- vcount(g)

    ppr_mat <- sapply(seq_len(n), function(i) {
    pers <- rep(0, n)
    pers[i] <- 1
    page_rank(g, personalized = pers, damping=0.7)$vector
    })


    hc_topkgraphs = hclust(as.dist(res$DIST), method="ward.D2")
    cl_topkgraphs = cutree(hc_topkgraphs, num_communities)

    hc_jaccard = hclust(as.dist(1-jaccard_sim), method="ward.D2")
    cl_jaccard = cutree(hc_jaccard, num_communities)

    hc_dice = hclust(as.dist(1-dice_sim), method="ward.D2")
    cl_dice = cutree(hc_dice, num_communities)

    hc_laplacian = hclust(dist(scale(emb)), method="ward.D2")
    cl_laplacian = cutree(hc_laplacian, num_communities)

    hc_node2vec = hclust(dist(scale(node2vec_emb)), method="ward.D2")
    cl_node2vec = cutree(hc_node2vec, num_communities)

    hc_ppr = hclust(as.dist(1-ppr_mat), method="ward.D2")
    cl_ppr = cutree(hc_ppr, num_communities)


    # ----------------------------------------------- #
    #hc_node2vec = hclust(dist(emb), method="ward.D")
    #cl_node2vec = cutree(hc_node2vec, num_communities)
    #ids = as.numeric(names(cl_node2vec))
    ids = match(1:n_nodes, node_order)
    cl_node2vec = cl_node2vec[ids]
    #print(cl_node2vec)

    


    library(aricode)
    #membership_gt <- rep(1:num_communities, each = nodes_per_community)
    membership_gt = V(g)$community  #rep(1:length(sizes), times = sizes)

    ari_topkgraphs = ARI(membership_gt, cl_topkgraphs)
    ari_jaccard = ARI(membership_gt, cl_jaccard)
    ari_dice = ARI(membership_gt, cl_dice)
    ari_laplacian = ARI(membership_gt, cl_laplacian)
    ari_ppr = ARI(membership_gt, cl_ppr)
    

    if(any(is.na(cl_node2vec))){
      #ari_node2vec  = NaN  
      na_id = which(is.na(cl_node2vec))
      cl_node2vec[na_id] = length(cl_node2vec) + na_id
      #print(cl_node2vec)
      ari_node2vec = ARI(membership_gt, cl_node2vec) 
    }else{
      ari_node2vec = ARI(membership_gt, cl_node2vec)
    }
    
    RES[xx,1] = ari_topkgraphs
    RES[xx,2] = ari_jaccard
    RES[xx,3] = ari_dice
    RES[xx,4] = ari_laplacian
    RES[xx,5] = ari_ppr
    RES[xx,6] = ari_node2vec
    

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

