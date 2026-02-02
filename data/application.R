# Citeseer
#edge_mat = read.table("Citeseer/citeseer_edge_matrix.csv", sep=",")
#node_labels = read.table("Citeseer/citeseer_node_labels.csv")
#node_labels = unlist(node_labels)

library(igraph)

# Paths
content_path <- "cora_data/cora/cora.content"
cites_path   <- "cora_data/cora/cora.cites"

# ---------------------------
# 1. Load nodes & labels
# ---------------------------
nodes <- read.table(content_path, stringsAsFactors = FALSE, sep = "\t")
node_ids <- nodes[, 1]               # paper IDs
labels_true <- as.factor(nodes[, ncol(nodes)])  # last column = class label

# Map node IDs -> consecutive integers (igraph vertex IDs)
id_map <- setNames(seq_along(node_ids), node_ids)

# ---------------------------
# 2. Load edges & remap IDs
# ---------------------------
edges_raw <- read.table(cites_path, stringsAsFactors = FALSE, sep = "")
colnames(edges_raw) <- c("from", "to")

edges_mapped <- data.frame(
  from = id_map[as.character(edges_raw$from)],
  to   = id_map[as.character(edges_raw$to)]
)

# Remove edges with missing nodes (external nodes)
edges_mapped <- na.omit(edges_mapped)

# ---------------------------
# 3. Build undirected igraph
# ---------------------------
g <- graph_from_edgelist(as.matrix(edges_mapped), directed = FALSE)

# ---------------------------
# 4. Attach labels to graph vertices
# ---------------------------
V(g)$label <- labels_true  # now labels align with vertices 1:n

# ---------------------------
# 5. Report
# ---------------------------
cat("Nodes:", vcount(g), "\nEdges:", ecount(g), "\nClasses:", 
    length(unique(labels_true)), "\n")



node_labels = labels_true
edge_mat = edges_mapped


n_runs = 30

RES_list = list()
ARI_res = matrix(NaN, n_runs, 2)
colnames(ARI_res) = c("TopKGraphs","Node2Vec")

NMI_res = matrix(NaN, n_runs, 2)
colnames(NMI_res) = c("TopKGraphs","Node2Vec")


for(yy in 1:n_runs){

    library(igraph)
    g = graph_from_edgelist(as.matrix(edge_mat), directed = FALSE)

    library(TopKLists)
    library(igraph)
    library(fastcluster) 
    library(cluster) 

    source("/home/bastian/GitHub/TopKGraphs/R/calc_SIL.R")
    source("/home/bastian/GitHub/TopKGraphs/R/calc_BINARY.R")
    source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs.R")
    source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs_walk.R")
    source("/home/bastian/GitHub/TopKGraphs/simulations/call_kNN.R")

    # Assign community membership
    V(g)$community = node_labels

    # Identify all connected components
    #components <- components(g)

    # Extract the largest component
    #largest <- induced_subgraph(g, which(components$membership == which.max(components$csize)))

    n_sub <- 50

    #print("Where is the error?")
    nodes = 1
    while((length(unique(nodes))) < n_sub){
            nodes <- random_walk(g, start = sample(V(g), 1),
                     steps = n_sub * 100,
                     mode = "all")
    }
    #print("Äh")
    #print(nodes)
    #print(unique(nodes))
    
    nodes <- unique(nodes)[1:n_sub]

    #print("YO")
    #print(nodes)

    # Induced subgraph
    largest <- induced_subgraph(g, nodes)

    V(largest)$community = droplevels(V(largest)$community)


    print(table(V(largest)$community))
    

    #plot(sub_g)


    cat("Nodes:", vcount(largest), "\nEdges:", ecount(largest), "\nClasses:", 
    length(unique(labels_true)), "\n")


    ## CALL TopKGraphs
    ###################################
    res = topkgraphs(list(largest), walk_depth=30, n_iter=50, n_cores=NaN)

    # clustering
    hc = hclust(as.dist(res$DIST), method="ward.D2")
    cl = cutree(hc, length(unique(V(largest)$community)))

    library(aricode)
    ari_topkgraphs = ARI(cl, V(largest)$community)
    nmi_topkgraphs = NMI(cl, V(largest)$community)

    ARI_res[yy,1] =  ari_topkgraphs
    NMI_res[yy,1] =  nmi_topkgraphs

    # classification
   # knn_topkgraphs = call_kNN_dist(res$DIST, V(largest)$community, k = 5)


    #library(caret)
    #all_levels <- sort(unique(c(knn_topkgraphs[[1]],knn_topkgraphs[[2]])))
    #cm_topkgraphs = confusionMatrix(factor(knn_topkgraphs[[1]], levels=all_levels), 
    #                                factor(knn_topkgraphs[[2]], levels=all_levels))


    #acc_topkgraphs = cm_topkgraphs$overall["Accuracy"]
    #bacc_topkgraphs = mean(cm_topkgraphs$byClass[, "Balanced Accuracy"],
    #                                na.rm = TRUE)

    ## CALL Node2Vec
    ###################################
    source("/home/bastian/GitHub/TopKGraphs/simulations/call_node2vec.R")
    node2vec_emb = call_node2vec(largest, walk_length=30, num_walks=50)

    #print(dim(node2vec_emb))
    node_order = round(as.numeric(rownames(node2vec_emb)))

    #print(node_order)
    node2vec_emb = matrix(as.numeric(unlist(node2vec_emb)), 
                                nrow=nrow(node2vec_emb), 
                                ncol=ncol(node2vec_emb))

    ids = match(1:length(V(largest)), node_order)

    # clustering
    hc_node2vec = hclust(dist(scale(node2vec_emb[ids,-1])), method="ward.D2")
    cl_node2vec = cutree(hc_node2vec,  length(unique(V(largest)$community)))

    # ----------------------------------------------- #
    #hc_node2vec = hclust(dist(emb), method="ward.D")
    #cl_node2vec = cutree(hc_node2vec, num_communities)
    #ids = as.numeric(names(cl_node2vec))

    #cl_node2vec = cl_node2vec[ids]

    ari_node2vec = ARI(cl_node2vec, V(largest)$community)
    nmi_node2vec = NMI(cl_node2vec, V(largest)$community)

    ARI_res[yy,2] =  ari_node2vec
    NMI_res[yy,2] =  nmi_node2vec


    # classification 
    #knn_node2vec = call_kNN_dist(as.matrix(dist(scale(node2vec_emb[ids,-1]))),
    #                                    V(largest)$community, k = 5)

    #all_levels <- sort(unique(c(knn_node2vec[[1]],knn_node2vec[[2]])))
    #    cm_node2vec = confusionMatrix(factor(knn_node2vec[[1]], levels=all_levels), 
    #                                    factor(knn_node2vec[[2]], levels=all_levels))
                                        
    #acc_node2vec = cm_node2vec$overall["Accuracy"]
    #bacc_node2vec <- mean(cm_node2vec$byClass[, "Balanced Accuracy"], na.rm = TRUE)

    ### CHECK the performance with varying k in kNN.

   # K = c(3,5,10,20,30,50)

    #RES = matrix(NaN,2,length(K))
    #rownames(RES) = c("TopKGraphs","Node2Vec")
    #colnames(RES) = K

   # for(xx in 1:length(K)){

        # Call TopKgraphs
        #knn_topkgraphs = call_kNN_dist(res$DIST, V(largest)$community, k = K[xx])
        
        #all_levels <- sort(unique(c(knn_topkgraphs[[1]],knn_topkgraphs[[2]])))
        #cm_topkgraphs = confusionMatrix(factor(knn_topkgraphs[[1]], levels=all_levels), 
         #                               factor(knn_topkgraphs[[2]], levels=all_levels))

        #bacc_topkgraphs = mean(cm_topkgraphs$byClass[, "Balanced Accuracy"],
        #                                na.rm = TRUE)

        #RES[1, xx] = bacc_topkgraphs

        # Call Node2Vec
        #knn_node2vec = call_kNN_dist(as.matrix(dist(scale(node2vec_emb[ids,-1]))),
        #                                V(largest)$community, k = K[xx])

        #all_levels <- sort(unique(c(knn_node2vec[[1]],knn_node2vec[[2]])))
        #cm_node2vec = confusionMatrix(factor(knn_node2vec[[1]], levels=all_levels), 
        #                                    factor(knn_node2vec[[2]], levels=all_levels))
                                            
        #bacc_node2vec <- mean(cm_node2vec$byClass[, "Balanced Accuracy"], na.rm = TRUE)

        #RES[2, xx] = bacc_node2vec

    #print(RES)

    #}

#RES_list[[yy]] = RES

#print(RES_list)
print("ARI")
print(ARI_res)
print("NMI")
print(NMI_res)
}
### Plots