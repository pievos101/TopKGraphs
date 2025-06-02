library(TopKLists)
library(igraph)


topkgraphs <- function(views, walk_depth=20, n_iter=20){

 RES = list()
 
 if(!(is(views[[1]])[1]=="igraph")){
    n_nodes = nrow(views[[1]])
 }else{
    n_nodes = length(V(views[[1]]))
 }

 cat("The graph has", n_nodes, "Nodes. \n")
 for(xx in 1:n_nodes){

    RES[[xx]] = topkgraphs_walk(views, 
                            start_node=xx, 
                            walk_depth=walk_depth, 
                            n_iter=n_iter)
 }

#return(RES)

# Convert to distance matrix
DIST = matrix(NaN, length(RES), length(RES))

for (xx in 1:length(RES)){

    DIST[xx, as.numeric(RES[[xx]]$TopK[,1])] = RES[[xx]]$Scores[,1]

}

DIST = DIST/max(DIST, na.rm=TRUE)
DIST[is.na(DIST)] = 1
diag(DIST) = 0

return(list(RES=RES, DIST=DIST))

}



topkgraphs_walk <- function(views, start_node=1, walk_depth=20, n_iter=20){

WALKSALL = list()

for(ii in 1:length(views)){

    adjmatrix = views[[ii]]

    if(!(is(adjmatrix)[1]=="igraph")){
        # Convert to graph object
        graph = graph_from_adjacency_matrix(
        adjmatrix,
        mode = "undirected",
        weighted = NULL,
        diag = FALSE,
        add.colnames = NULL,
        add.rownames = NA
        )
    }else{
        graph = adjmatrix
    }

    WALKS = matrix(NaN, walk_depth + 1, n_iter)

        for(xx in 1:n_iter){

            # Generate Random walks
            list = random_walk(
            graph,
            start_node,
            walk_depth,
            mode = "all", #c("out", "in", "all", "total"),
            stuck = c("return")
            )

            list = as.numeric(list)
            WALKS[,xx] = list
        }

    WALKSALL[[ii]] = WALKS

}# End of iterations over graphs


# Concatenate the walks from the multi-view graphs
WALKS = Reduce('cbind', WALKSALL) 

#print(WALKS)

res1 = apply(WALKS, 2, table, simplify=FALSE)

#print(res1)

res2 = sapply(res1, sort, decreasing=TRUE, simplify=FALSE)

#print(res2)
#@FIXME ?? - this should be fixed now
#n_nodes = max(as.numeric(unlist(sapply(res2, names))))
n_nodes = length(unique(as.numeric(unlist(sapply(res2, names)))))
#print(n_nodes)
#print(ncol(WALKS))

rankMatrix = matrix(NaN, n_nodes, ncol(WALKS))

# Fill the rank matrix
for (xx in 1:ncol(WALKS)){
    
    rankMatrix[1:length(res2[[xx]]),xx] = names(res2[[xx]])

}

#print(rankMatrix)

# check missings
values = unique(rankMatrix)
na_ids = which(rankMatrix=="NaN", arr.ind=TRUE)
na_cols = na_ids[,2]

for(xx in 1:length(na_cols)){

   repl_val =  setdiff(values, rankMatrix[, na_cols[xx]])
   #print(rankMatrix[, na_cols[xx]])
   #print(repl_val)
   if(length(repl_val)==0){next}
   repl_ids = which(rankMatrix[,na_cols[xx]]=="NaN")
   if(length(repl_ids)==0){next}
   #print(repl_val)
   #print("start")
   #print(na_cols[xx])
   #print(rankMatrix)
   rankMatrix[repl_ids, na_cols[xx]] = sample(repl_val, length(repl_val), 
                                                replace=FALSE) 
   #xx = xx + length(repl_ids) - 1
   #print("done")

}

### check for remaining NaNs

#print(rankMatrix)

na_ids = which(rankMatrix=="NaN", arr.ind=TRUE)
na_rows = unique(na_ids[,1])

#print(na_rows)

if(length(na_rows)!=0){
    rankMatrix = rankMatrix[-na_rows,]
}

#print(rankMatrix)

# Check the correctness of that! @FIXME
# BORDA ---------------------------------------- #
#print("BORDA")
rankMatrix2 <- apply(rankMatrix,2, function(x) match(1:length(x),x))

#print(rankMatrix2)

rankMatrix2 <- matrix(as.character(rankMatrix2), dim(rankMatrix2)[1],dim(rankMatrix2)[2])

#print(rankMatrix2)

rankMatrix2 = rankMatrix

IN          <- lapply(seq_len(ncol(rankMatrix2)), function(i) rankMatrix2[,i])

#print(IN)

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

return(borda.res)

}