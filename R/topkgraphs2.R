library(igraph)
library(foreach)
library(doParallel)
library(Rcpp)

# Source your C++ biased walk function
Rcpp::sourceCpp("/home/bpfeif/GitHub/TopKGraphs/src/walk_with_jaccard_degree_safe.cpp")

topkgraphs2 <- function(views, 
                        walk_depth = 50, 
                        n_iter = 50, 
                        n_cores = NA, 
                        alpha = 0,     # restart probability
                        lambda = 0.3     # decay factor per first-hit step
                        ){
  
  # -------------------------------
  # Internal function: generate walks and compute per-walk first-hit decay scores
  # -------------------------------
  topkgraphs_walk2 <- function(views, start_node, walk_depth, n_iter, alpha){
    
    n_nodes <- if(!(is(views[[1]])[1]=="igraph")) nrow(views[[1]]) else length(V(views[[1]]))
    decay_scores <- matrix(0, nrow = n_nodes, ncol = n_iter * length(views))
    
    col_counter <- 1
    
    for(ii in seq_along(views)){
      adjmatrix <- views[[ii]]
      
      if(!(is(adjmatrix)[1] == "igraph")){
        graph <- graph_from_adjacency_matrix(adjmatrix, mode = "undirected")
      } else {
        graph <- adjmatrix
      }
      
      adj_list <- lapply(V(graph), function(v) as.integer(neighbors(graph, v)))
      
      for(iter in 1:n_iter){
        walk <- walk_with_jaccard_degree_safe_cpp(adj_list, start_node, walk_depth, alpha)
        
        # First-hit step per node
        first_hit <- rep(walk_depth + 1, n_nodes)
        for(pos in seq_along(walk)){
          node <- walk[pos]
          if(pos < first_hit[node]){
            first_hit[node] <- pos
          }
        }
        
        # Convert to decay score
        decay_scores[, col_counter] <- exp(-lambda * first_hit)
        col_counter <- col_counter + 1
      }
    }
    
    return(decay_scores)
  }
  
  # -------------------------------
  # Main function: compute per-node distances
  # -------------------------------
  n_nodes <- if(!(is(views[[1]])[1]=="igraph")) nrow(views[[1]]) else length(V(views[[1]]))
  cat("Graph has", n_nodes, "nodes.\n")
  
  RES <- vector("list", n_nodes)
  
  node_indices <- 1:n_nodes
  if(!is.na(n_cores)){
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    RES <- foreach(xx = node_indices, .combine = c, .packages = c("igraph")) %dopar% {
      list(topkgraphs_walk2(views, start_node = xx, walk_depth = walk_depth, n_iter = n_iter, alpha = alpha))
    }
    
    stopCluster(cl)
  } else {
    for(xx in node_indices){
      RES[[xx]] <- topkgraphs_walk2(views, start_node = xx, walk_depth = walk_depth, n_iter = n_iter, alpha = alpha)
    }
  }
  
  # -------------------------------
  # Aggregate per-node scores across walks
  # -------------------------------
  DIST <- matrix(0, n_nodes, n_nodes)
  
  for(xx in seq_len(n_nodes)){
    decay_scores <- RES[[xx]]  # n_nodes x total_walks
    # Aggregate per node (mean across walks)
    mean_scores <- rowMeans(decay_scores)
    DIST[xx, ] <- mean_scores
  }
  
  # Convert to distance
  DIST <- 1 - DIST / max(DIST)
  diag(DIST) <- 0
  DIST_raw <- (DIST + t(DIST)) / 2  # symmetrize
  
  # Optional embedding
  coords <- cmdscale(as.dist(DIST_raw), k = 10)
  DIST_embedded <- as.matrix(dist(coords))
  
  return(list(RES = RES, DIST = DIST_embedded, DIST_raw = DIST_raw))
}
