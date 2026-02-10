#' Jaccard-Anchored Random Walk Ranking
#'
#' Performs multiple biased random walks from a start node and aggregates
#' the resulting partial rankings.
#'
#' @param views List of igraph objects or adjacency matrices
#' @param start_node Integer start node
#' @param walk_depth Length of each walk
#' @param n_iter Number of walks
#'
#' @return A list containing aggregated rankings and scores
#'
#' @export

topkgraphs_walk <- function(views, 
                            start_node=1, 
                            walk_depth=20, 
                            n_iter=20 
                            ){



do.BORDA = TRUE
do.TopKSignal = FALSE
do.RRA = FALSE
alpha=0

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


    WALKS <- walk_with_jaccard_degree_safe_cpp(adj_list,
                                           start_node = start_node,
                                           walk_depth = walk_depth,
                                           n_iter = n_iter)     

    WALKSALL[[ii]] = WALKS

}# End of iterations over graphs


# Concatenate the walks from the multi-view graphs
WALKS = Reduce('cbind', WALKSALL) 

#print(WALKS)

#res1 = apply(WALKS, 2, table, simplify=FALSE)

#print(res1)

#res2 = sapply(res1, sort, decreasing=TRUE, simplify=FALSE)

res3 = apply(WALKS, 2, unique, simplify=FALSE)


#print(res2)
#n_nodes1 = length(unique(as.numeric(unlist(sapply(res2, names)))))
n_nodes2 = length(unique(unlist(res3)))
#print(unique(unlist(res2)))
#n_nodes = unique(c(n_nodes1, n_nodes2))
#print(n_nodes)
#print(ncol(WALKS))

#rankMatrix1 = matrix(NaN, n_nodes1, ncol(WALKS))
#print(rankMatrix)
rankMatrix2 = matrix(NaN, n_nodes2, ncol(WALKS))

# Fill the rank matrix
for (xx in 1:ncol(WALKS)){
    
    #rankMatrix1[1:length(res2[[xx]]),xx] = names(res2[[xx]])
    rankMatrix2[1:length(res3[[xx]]),xx] = res3[[xx]]

}

rankMatrix = rankMatrix2 # cbind(rankMatrix1, rankMatrix2)

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
   #print(repl_val)

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

   if(length(repl_val)==1){
   repl_val_perm = repl_val 
   }else{
   repl_val_perm = sample(repl_val, length(repl_val), replace=FALSE) 
   }

   #print(repl_val_perm)
   
   #rankMatrix[repl_ids, na_cols[xx]] = repl_val_perm
   
   #xx = xx + length(repl_ids) - 1
   #print("done")
   #print("Updated Matrix")
  # print(rankMatrix)

}

### check for remaining NaNs

#print(rankMatrix)

#print(rankMatrix)
#print(dim(rankMatrix))
#print(length(unique(as.vector(rankMatrix))))
#print(unique(as.vector(rankMatrix)))


#na_ids = which(rankMatrix=="NaN", arr.ind=TRUE)
#na_rows = unique(na_ids[,1])

#print(na_rows)

#if(length(na_rows)!=0){
#    rankMatrix = rankMatrix[-na_rows,]
    #print("There are remaining NaN's")
#}

#print(rankMatrix)

# Check the correctness of that! @FIXME
# BORDA ---------------------------------------- #
#print("BORDA")
#rankMatrix2 <- apply(rankMatrix,2, function(x) match(1:length(x),x))

#print(rankMatrix2)

#rankMatrix2 <- matrix(as.character(rankMatrix2), dim(rankMatrix2)[1],dim(rankMatrix2)[2])

#print(rankMatrix2)

rankMatrix2 = rankMatrix
borda.res = NaN

if(do.BORDA){

IN          <- lapply(seq_len(ncol(rankMatrix2)), function(i) rankMatrix2[,i])
#print(IN)
#print(IN)

IN = lapply(IN, function(x){as.numeric(x[!is.na(as.numeric(x))])})
#print(IN)
require(TopKLists)
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

}


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
latent_signal = NaN
if(do.TopKSignal){


  node_names = sort(unique(as.vector(rankMatrix2)))
  IN = apply(rankMatrix2, 2, function(x){match(node_names,x)})
  rownames(IN) = node_names
  colnames(IN) = paste("r",1:ncol(rankMatrix2), sep="")


  require(TopKSignal)
  require(gurobi)

  estimatedSignal <- estimateTheta(R.input = IN, 
          num.boot = 50, b = 0.1, solver = "gurobi", 
          type = "restrictedQuadratic", bootstrap.type = "classic.bootstrap",
          nCore = 5)

  latent_signal = estimatedSignal$estimation

}

agg.res = list()
if(do.RRA){
  require(RobustRankAggreg)
  node_names = sort(unique(as.vector(rankMatrix2)))
  IN = apply(rankMatrix2, 2, function(x){match(node_names,x)})
  rownames(IN) = node_names
  colnames(IN) = paste("r",1:ncol(rankMatrix2), sep="")
  rankMatrix = IN
  rankMatrix2 <- apply(rankMatrix,2, function(x) match(1:length(x),x))
  rankMatrix2 <- matrix(as.character(rankMatrix2), dim(rankMatrix2)[1],dim(rankMatrix2)[2])
  IN          <- lapply(seq_len(ncol(rankMatrix2)), function(i) rankMatrix2[,i])
  agg.res     <- aggregateRanks(IN)

}

return(list(BORDA=borda.res, RM=rankMatrix2, TS=latent_signal, RRA=agg.res))

}