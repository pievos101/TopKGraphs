# ======================================================
# 0. Libraries
# ======================================================
library(igraph)
library(ggplot2)
library(reshape2)
library(aricode)   # ARI, NMI, AMI
library(caret)     # confusionMatrix

# ======================================================
# 1. Load TopKGraphs functions
# ======================================================
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bastian/GitHub/TopKGraphs/simulations/call_node2vec.R")
source("/home/bastian/GitHub/TopKGraphs/simulations/call_kNN.R")

# ======================================================
# 2. External metrics helper
# ======================================================
compute_external_metrics <- function(labels_true, labels_pred) {
  c(
    ARI = ARI(labels_true, labels_pred),
    NMI = NMI(labels_true, labels_pred),
    AMI = AMI(labels_true, labels_pred)
  )
}

# ======================================================
# 3. Load dataset: CORA
# ======================================================
DATASET <- "CORA"

if(DATASET=="CORA"){

  # Paths
  content_path <- "cora_data/cora/cora.content"
  cites_path   <- "cora_data/cora/cora.cites"

  # Load nodes & labels
  nodes <- read.table(content_path, stringsAsFactors = FALSE, sep = "\t")
  node_ids <- nodes[,1]
  labels_true <- as.factor(nodes[,ncol(nodes)])

  # Map node IDs -> consecutive integers
  id_map <- setNames(seq_along(node_ids), node_ids)

  # Load edges
  edges_raw <- read.table(cites_path, stringsAsFactors = FALSE, sep = "")
  colnames(edges_raw) <- c("from","to")
  edges_mapped <- data.frame(
    from = id_map[as.character(edges_raw$from)],
    to   = id_map[as.character(edges_raw$to)]
  )
  edges_mapped <- na.omit(edges_mapped)

  # Build graph
  g <- graph_from_edgelist(as.matrix(edges_mapped), directed = FALSE)
  V(g)$label <- labels_true

  cat("Nodes:", vcount(g), "Edges:", ecount(g), "Classes:", length(unique(labels_true)), "\n")

  node_labels <- labels_true
  edge_mat <- edges_mapped
}

# ======================================================
# 4. Experiment setup
# ======================================================
n_runs <- 50
methods <- c("TopKGraphs","Node2Vec","Jaccard","Dice","PageRank")
metrics <- c("ARI","NMI","AMI")
K_vals <- c(5,7,10)
n_nodes_sub <- 200

# Arrays to store results
CLUSTER_array <- array(NA, dim=c(n_runs, length(methods), length(metrics)),
                       dimnames=list(NULL, methods, metrics))
KNN_array     <- array(NA, dim=c(n_runs, length(methods), length(K_vals)),
                       dimnames=list(NULL, methods, paste0("K=",K_vals)))

# ======================================================
# 5. Main loop
# ======================================================
for(yy in 1:n_runs){

  cat("Run:", yy, "\n")

  # -----------------------
  # 5a. Subgraph sampling
  # -----------------------
  nodes_sample <- 1
  while(length(unique(nodes_sample)) < n_nodes_sub){
    nodes_sample <- random_walk(g, start=sample(V(g),1),
                                steps=n_nodes_sub*100, mode="all")
  }
  nodes_sample <- unique(nodes_sample)[1:n_nodes_sub]
  sub_g <- induced_subgraph(g, nodes_sample)
  V(sub_g)$community <- droplevels(V(sub_g)$label)
  all_levels <- levels(V(sub_g)$community)  # for kNN factors

  # -----------------------
  # 5b. TopKGraphs
  # -----------------------
  res_tkg <- topkgraphs(list(sub_g), walk_depth=100, n_iter=50, n_cores=NaN)
  hc_tkg <- hclust(as.dist(res_tkg$DIST), method="ward.D2")
  cl_tkg <- cutree(hc_tkg, length(unique(V(sub_g)$community)))

  # -----------------------
  # 5c. Node2Vec
  # -----------------------
  node2vec_emb <- call_node2vec(sub_g, walk_length=100, num_walks=50)
  node_order <- round(as.numeric(rownames(node2vec_emb)))
  emb <- matrix(as.numeric(unlist(node2vec_emb)), nrow=nrow(node2vec_emb))
  ids <- match(1:vcount(sub_g), node_order)
  hc_n2v <- hclust(dist(scale(emb[ids,-1])), method="ward.D2")
  cl_n2v <- cutree(hc_n2v, length(unique(V(sub_g)$community)))

  # -----------------------
  # 5d. Jaccard & Dice
  # -----------------------
  jacc_sim <- similarity(sub_g, method="jaccard", mode="all")
  dice_sim <- similarity(sub_g, method="dice", mode="all")
  cl_jacc <- cutree(hclust(as.dist(1-jacc_sim), method="ward.D2"),
                    length(unique(V(sub_g)$community)))
  cl_dice <- cutree(hclust(as.dist(1-dice_sim), method="ward.D2"),
                    length(unique(V(sub_g)$community)))

  # -----------------------
  # 5e. PageRank
  # -----------------------
  n <- vcount(sub_g)
  ppr_mat <- sapply(1:n, function(i){
    pers <- rep(0,n); pers[i] <- 1
    page_rank(sub_g, personalized=pers, damping=0.7)$vector
  })
  cl_ppr <- cutree(hclust(as.dist(1-ppr_mat), method="ward.D2"),
                   length(unique(V(sub_g)$community)))

  # -----------------------
  # 6. Store clustering metrics
  # -----------------------
  CLUSTER_array[yy,"TopKGraphs",] <- compute_external_metrics(V(sub_g)$community, cl_tkg)
  CLUSTER_array[yy,"Node2Vec",] <- compute_external_metrics(V(sub_g)$community, cl_n2v)
  CLUSTER_array[yy,"Jaccard",] <- compute_external_metrics(V(sub_g)$community, cl_jacc)
  CLUSTER_array[yy,"Dice",] <- compute_external_metrics(V(sub_g)$community, cl_dice)
  CLUSTER_array[yy,"PageRank",] <- compute_external_metrics(V(sub_g)$community, cl_ppr)

  print(CLUSTER_array)  

  # -----------------------
  # 7. Compute kNN Balanced Accuracy for all K & all methods
  # -----------------------
  safe_balacc <- function(knn_result){
  # Extract from the list
  y_pred <- knn_result$y_pred
  y_test <- knn_result$y_test

  # Check if at least 2 classes
  if(length(unique(y_pred)) < 2 || length(unique(y_test)) < 2) {
    return(NA)
  }

  # Make factors with matching levels
  levels_all <- sort(unique(c(y_pred, y_test)))
  pred <- factor(y_pred, levels = levels_all)
  ref  <- factor(y_test, levels = levels_all)

  cm <- confusionMatrix(pred, ref)
  
  # Compute Balanced Accuracy
  if(is.null(dim(cm$byClass)[1])){
    bacc <- cm$byClass["Balanced Accuracy"]
  } else {
    bacc <- mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE)
  }
  return(bacc)
  }
  # -------------------------------------------- #

  for(xx in seq_along(K_vals)){
    k_nn <- K_vals[xx]

    # TopKGraphs
    knn_tkg <- call_kNN_dist(res_tkg$DIST, V(sub_g)$community, k=k_nn)
    KNN_array[yy,"TopKGraphs",xx] <- safe_balacc(knn_tkg)

    # Node2Vec
    knn_n2v <- call_kNN_dist(as.matrix(dist(scale(emb[ids,-1]))),
                             V(sub_g)$community, k=k_nn)
    KNN_array[yy,"Node2Vec",xx] <- safe_balacc(knn_n2v)

    # Jaccard
    knn_jacc <- call_kNN_dist(as.matrix(dist(1-jacc_sim)),
                              V(sub_g)$community, k=k_nn)
    KNN_array[yy,"Jaccard",xx] <- safe_balacc(knn_jacc)

    # Dice
    knn_dice <- call_kNN_dist(as.matrix(dist(1-dice_sim)),
                              V(sub_g)$community, k=k_nn)
    KNN_array[yy,"Dice",xx] <- safe_balacc(knn_dice)

    # PageRank
    knn_ppr <- call_kNN_dist(as.matrix(dist(1-ppr_mat)),
                             V(sub_g)$community, k=k_nn)
    KNN_array[yy,"PageRank",xx] <- safe_balacc(knn_ppr)
  }

print(KNN_array)

} # end runs

library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)  # for viewport

# ======================================================
# Prepare kNN plot
# ======================================================
df_knn <- melt(KNN_array)
colnames(df_knn) <- c("Run","Method","K","BalancedAccuracy")
df_knn$K <- as.factor(df_knn$K)

p_knn <- ggplot(df_knn, aes(x=K, y=BalancedAccuracy, fill=Method)) +
  geom_boxplot(position=position_dodge(width=0.75), width=0.6, outlier.alpha=0.6) +
  labs(title="", x="k", y="Balanced Accuracy") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(face="bold"),
        legend.position="top")  # keep legend here

# ======================================================
# Prepare clustering metrics plot
# ======================================================
df_cluster <- melt(CLUSTER_array)
colnames(df_cluster) <- c("Run","Method","Metric","Score")

p_cluster <- ggplot(df_cluster, aes(x=Method, y=Score, fill=Method)) +
  geom_bar(stat="summary", fun="mean", position="dodge") +
  geom_errorbar(stat="summary", fun.data="mean_se",
                width=0.2, position=position_dodge(0.9)) +
  facet_wrap(~Metric, scales="free_y") +
  theme_minimal(base_size=14) +
  theme(
    plot.title = element_text(face="bold"),
    legend.position="none",       # remove legend here
    axis.text.x = element_blank() # remove x-axis labels
  ) +
  labs(title="", x="Method", y="Score")

# ======================================================
# Extract shared legend from p_knn
# ======================================================
get_legend <- function(a_plot){
  tmp <- ggplotGrob(a_plot)
  leg_index <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg_index]]
  return(legend)
}

shared_legend <- get_legend(p_knn)

# Remove legend from p_knn
p_knn <- p_knn + theme(legend.position="none")

# ======================================================
# Arrange side by side with shared legend on top
# ======================================================
grid.arrange(
  shared_legend,
  arrangeGrob(p_knn, p_cluster, ncol=2),
  heights=c(0.15, 0.85)  # adjust space for legend
)