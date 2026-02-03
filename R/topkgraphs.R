library(TopKLists)
#library(TopKSignal)
library(igraph)
#source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs_walk.R")

topkgraphs <- function(views, walk_depth=50, n_iter=50, 
                                 n_cores=NaN, agg_method="mean",
                                 do.BORDA = TRUE, 
                                 do.TopKSignal=FALSE,
                                 do.RRA=FALSE){

 RES = list()
 
 if(!(is(views[[1]])[1]=="igraph")){
    n_nodes = nrow(views[[1]])
 }else{
    n_nodes = length(V(views[[1]]))
 }

 cat("The graph has", n_nodes, "Nodes. \n")
 
 if(is.na(n_cores)){ # no parallel computation
    for(xx in 1:n_nodes){

        RES[[xx]] = topkgraphs_walk(views, 
                                start_node=xx, 
                                walk_depth=walk_depth, 
                                n_iter=n_iter,
                                do.BORDA=do.BORDA, 
                                do.TopKSignal=do.TopKSignal,
                                do.RRA=do.RRA)
    }
 }else{ # parallel computation
 
    require(foreach)
    require(doParallel)

    ncores <- parallel::detectCores() - 3
    cl <- makeCluster(ncores)

    # Register the cluster
    registerDoParallel(cl)

    require(Rcpp)
    Rcpp::sourceCpp("/home/bastian/GitHub/TopKGraphs/src/walk_with_jaccard_degree_safe.cpp")

    RES = foreach(xx = 1:n_nodes, .combine = c,
                    .export = "topkgraphs_walk",
                    .packages = c("igraph","TopKLists") ) %dopar% {
                     list(topkgraphs_walk(views, 
                                start_node=xx, 
                                walk_depth=walk_depth, 
                                n_iter=n_iter,
                                do.BORDA=do.BORDA, 
                                do.TopKSignal=do.TopKSignal,
                                do.RRA=do.RRA))
    }
 }

#return(RES)

#print(RES)

# Convert to distance matrix
DIST = matrix(NaN, length(RES), length(RES))
#print(RES)
for (xx in 1:length(RES)){


   if(do.BORDA){

      if(agg_method=="mean"){
         DIST[xx, as.numeric(RES[[xx]]$BORDA$TopK[,1])] = 
         RES[[xx]]$BORDA$Scores[,1]
      }
      if(agg_method=="median"){
         DIST[xx, as.numeric(RES[[xx]]$BORDA$TopK[,2])] = 
         RES[[xx]]$BORDA$Scores[,2]
      }
      if(agg_method=="geo.mean"){
         DIST[xx, as.numeric(RES[[xx]]$BORDA$TopK[,3])] = 
         RES[[xx]]$BORDA$Scores[,3]
      }
      if(agg_method=="l2norm"){
         DIST[xx, as.numeric(RES[[xx]]$BORDA$TopK[,4])] = 
         RES[[xx]]$BORDA$Scores[,4]
      }
      
      }

   if(do.TopKSignal){
         DIST[xx, as.numeric(RES[[xx]]$TS$id)] = 
         1 / (1 + exp(-RES[[xx]]$TS$signal.estimate))
   }

   if(do.RRA){
      DIST[xx, as.numeric(RES[[xx]]$RRA$Name)] = 
         1 / (1 + exp(-RES[[xx]]$RRA$Score))
   }


}

if(do.BORDA){
#DIST = DIST/max(DIST, na.rm=TRUE)
}

if(do.TopKSignal){
DIST[DIST<0] = 0
DIST = 1 - DIST/max(DIST, na.rm=TRUE)
}

if(do.RRA){
DIST = DIST/max(DIST, na.rm=TRUE)
}

DIST[is.na(DIST)] <- max(DIST, na.rm=TRUE)
diag(DIST) <- 0
DIST_raw <- (DIST + t(DIST)) / 2  # symmetrize

coords <- cmdscale(as.dist(DIST_raw), k=10)  # 10D embedding
DIST <- as.matrix(dist(coords))

return(list(RES=RES, DIST=DIST, DIST_raw=DIST_raw))

}

