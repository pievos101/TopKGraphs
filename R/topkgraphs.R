library(TopKLists)
library(igraph)


adjmatrix1 = omics_binary[[1]]
adjmatrix2 = omics_binary[[2]]

# Convert to graph object
graph1 = graph_from_adjacency_matrix(
  adjmatrix1,
  mode = "undirected",
  weighted = NULL,
  diag = FALSE,
  add.colnames = NULL,
  add.rownames = NA
)

# Convert to graph object
graph2 = graph_from_adjacency_matrix(
  adjmatrix2,
  mode = "undirected",
  weighted = NULL,
  diag = FALSE,
  add.colnames = NULL,
  add.rownames = NA
)


walk_depth = 10 
n_iter = 20 
start_node = 1

WALKS1 = matrix(NaN, walk_depth + 1, n_iter)

for(xx in 1:n_iter){

    # Generate Random walks
    list = random_walk(
    graph1,
    start_node,
    walk_depth,
    mode = "all", #c("out", "in", "all", "total"),
    stuck = c("return")
    )

    list = as.numeric(list)
    WALKS1[,xx] = list
}


WALKS2 = matrix(NaN, walk_depth + 1, n_iter)

for(xx in 1:n_iter){

    # Generate Random walks
    list = random_walk(
    graph2,
    start_node,
    walk_depth,
    mode = "all", #c("out", "in", "all", "total"),
    stuck = c("return")
    )

    list = as.numeric(list)
    WALKS2[,xx] = list
}

# Concatenate the walks from the multi-view graphs
WALKS = cbind(WALKS1, WALKS2)

res1 = apply(WALKS, 2, table)
res2 = sapply(res1, sort, decreasing=TRUE)
n_nodes = max(sapply(res2, length))

rankMatrix = matrix(NaN, n_nodes, ncol(WALKS))

# Fill the rank matrix
for (xx in 1:ncol(WALKS)){
    
    rankMatrix[1:length(res2[[xx]]),xx] = names(res2[[xx]])

}

# check missings
values = unique(rankMatrix)
na_ids = which(rankMatrix=="NaN", arr.ind=TRUE)

for(xx in 1:nrow(na_ids)){

   repl_val =  setdiff(values, rankMatrix[,na_ids[xx,2]])
   #print(repl_val)
   repl_ids = which(rankMatrix[,na_ids[xx,2]]=="NaN")
   rankMatrix[repl_ids, na_ids[xx,2]] = sample(repl_val, length(repl_val), 
                                                replace=FALSE) 

}

# Check the correctness of that! @FIXME
# BORDA ---------------------------------------- #
print("BORDA")
rankMatrix2 <- apply(rankMatrix,2, function(x) match(1:length(x),x))
rankMatrix2 <- matrix(as.character(rankMatrix2), dim(rankMatrix2)[1],dim(rankMatrix2)[2])
IN          <- lapply(seq_len(ncol(rankMatrix2)), function(i) rankMatrix2[,i])
borda.res   <- Borda(IN)
# l2norm
TKSrank_l2norm   <- borda.res$TopK$l2norm
TKSrank_l2norm   <- match(1:length(TKSrank_l2norm),as.numeric(TKSrank_l2norm))
# mean
TKSrank_mean     <- borda.res$TopK$mean
TKSrank_mean     <- match(1:length(TKSrank_mean),as.numeric(TKSrank_mean))
# median
TKSrank_median   <- borda.res$TopK$median
TKSrank_median   <- match(1:length(TKSrank_median),as.numeric(TKSrank_median))
# geo.mean
TKSrank_geomean   <- borda.res$TopK$geo.mean
TKSrank_geomean   <- match(1:length(TKSrank_geomean),as.numeric(TKSrank_geomean))
#print(length(TKSrank))
# ---------------------------------------------- #