#
walk_with_restart_fast <- function(adj_list, start_node, walk_depth, alpha=0.3) {
  walk <- numeric(walk_depth + 1)
  walk[1] <- start_node
  current <- start_node
  
  for (i in 2:(walk_depth + 1)) {
    if (runif(1) < alpha) {
      current <- start_node  # restart
    } else {
      neighbors <- adj_list[[current]]
      if (length(neighbors) == 0) break
      current <- sample(neighbors, 1)
    }
    walk[i] <- current
  }
  return(walk)
}
#
walk_with_restart <- function(graph, start_node, walk_depth, alpha=0.3){

  walk <- numeric(walk_depth + 1)
  walk[1] <- start_node
  current <- start_node
  
  for (i in 2:(walk_depth+1)){
    if(runif(1) < alpha){
      current <- start_node
    } else {
      neighbors <- neighbors(graph, current)
      if(length(neighbors)==0) break
      current <- sample(neighbors, 1)
    }
    walk[i] <- current
  }
  
  return(walk)
}


topkgraphs_walk <- function(views, 
                            start_node=1, 
                            walk_depth=20, 
                            n_iter=20, 
                            alpha=0.3){

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
    
    # Precompute adjacency list
    adj_list <- lapply(V(graph), function(v) as.integer(neighbors(graph, v)))

        #for(xx in 1:n_iter){

            # Generate Random walks
        #    list = walk_with_restart_fast(
        #      adj_list,
        #      start_node,
        #      walk_depth
        #    )
            
            #list = walk_with_restart(
            #graph,
            #start_node,
            #walk_depth,
            #mode = "all", #c("out", "in", "all", "total"),
            #stuck = c("return")
            #)

         #   list = as.numeric(list)
            #print(list)
            #print(dim(WALKS))
         #   WALKS[,xx] = list
        #}

    WALKS = replicate(n_iter, walk_with_restart_fast(adj_list, 
                      start_node=start_node, walk_depth= walk_depth, 
                      alpha=alpha))      

    WALKSALL[[ii]] = WALKS

}# End of iterations over graphs


# Concatenate the walks from the multi-view graphs
WALKS = Reduce('cbind', WALKSALL) 

#print(WALKS)

#res1 = apply(WALKS, 2, table, simplify=FALSE)

#print(res1)

#res2 = sapply(res1, sort, decreasing=TRUE, simplify=FALSE)


res2 = apply(WALKS, 2, unique, simplify=FALSE)


#print(res2)
#@FIXME ?? - this should be fixed now
#n_nodes = max(as.numeric(unlist(sapply(res2, names))))
#n_nodes1 = length(unique(as.numeric(unlist(sapply(res2, names)))))
n_nodes = length(unique(unlist(res2)))
#n_nodes = unique(c(n_nodes1, n_nodes2))
#print(n_nodes)
#print(ncol(WALKS))

rankMatrix = matrix(NaN, n_nodes, ncol(WALKS))
#rankMatrix2 = matrix(NaN, n_nodes2, ncol(WALKS))

# Fill the rank matrix
for (xx in 1:ncol(WALKS)){
    
    #rankMatrix[1:length(res2[[xx]]),xx] = names(res2[[xx]])
    rankMatrix[1:length(res2[[xx]]),xx] = res2[[xx]]

}

#rankMatrix = cbind(rankMatrix1, rankMatrix2)

#print(rankMatrix)

# check missings
#print("Walks")
#print(start_node)
#print(rankMatrix)
values = unique(as.vector(rankMatrix))
#print("Unique values")
#print(values)
na_ids = which(rankMatrix=="NaN", arr.ind=TRUE)
na_cols = unique(na_ids[,2])
#print(na_cols)
#stop("Check")

for(xx in 1:length(na_cols)){

   repl_val =  setdiff(values, rankMatrix[, na_cols[xx]])
   #print(rankMatrix[, na_cols[xx]])

   #print("Replacers")
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
   #print("Updated Matrix")
  # print(rankMatrix)

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
#rankMatrix2 <- apply(rankMatrix,2, function(x) match(1:length(x),x))

#print(rankMatrix2)

#rankMatrix2 <- matrix(as.character(rankMatrix2), dim(rankMatrix2)[1],dim(rankMatrix2)[2])

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

#rownames(rankMatrix2) = paste("p", 1:nrow(rankMatrix2), sep="")
#colnames(rankMatrix2) = paste("r", 1:ncol(rankMatrix2), sep="")

#estimatedSignal <- estimateTheta(R.input = rankMatrix2, 
#                    num.boot = 100, b = 0.1, 
#                    solver = "gurobi",
#                    #type = "restrictedLinear", 
#                    type = "restrictedQuadratic", 
#                    bootstrap.type = "classic.bootstrap",
#                    nCore = 3)



#TKSrank_topksignal = rank(-estimatedSignal$estimation$signal.estimate)

return(list(BORDA=borda.res, RM=rankMatrix2))

}