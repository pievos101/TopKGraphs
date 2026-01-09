# Citeseer
edge_mat = read.table("Citeseer/citeseer_edge_matrix.csv", sep=",")
node_labels = read.table("Citeseer/citeseer_node_labels.csv")
node_labels = unlist(node_labels)



library(igraph)
g = graph_from_edgelist(as.matrix(edge_mat), directed = FALSE)

library(TopKLists)
library(igraph)
library(fastcluster) 
library(cluster) 

source("/home/bastian/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_BINARY.R")
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs.R")

# Assign community membership
V(g)$community = node_labels

# Identify all connected components
components <- components(g)

# Extract the largest component
largest <- induced_subgraph(g, which(components$membership == which.max(components$csize)))


## CALL TopKGraphs
###################################
res = topkgraphs(list(largest), walk_depth=50, n_iter=30, n_cores=5)

hc = hclust(as.dist(res$DIST), method="ward.D2")

cl = cutree(hc, length(unique(V(largest)$community)))

library(aricode)

ari_topkgraphs = ARI(cl, V(largest)$community)

## CALL Node2Vec
###################################
source("/home/bastian/GitHub/TopKGraphs/simulations/call_node2vec.R")
node2vec_emb = call_node2vec(largest)

#print(dim(node2vec_emb))
node_order = round(as.numeric(rownames(node2vec_emb)))
#print(node_order)
node2vec_emb = as.matrix(node2vec_emb)

hc_node2vec = hclust(dist(scale(node2vec_emb)), method="ward.D2")
cl_node2vec = cutree(hc_node2vec,  length(unique(V(largest)$community)))

# ----------------------------------------------- #
#hc_node2vec = hclust(dist(emb), method="ward.D")
#cl_node2vec = cutree(hc_node2vec, num_communities)
#ids = as.numeric(names(cl_node2vec))
ids = match(1:length(V(largest)), node_order)
cl_node2vec = cl_node2vec[ids]

ari_node2vec = ARI(cl_node2vec, V(largest)$community)