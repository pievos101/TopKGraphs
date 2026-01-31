# ======================================================
# TopKGraphs kNN classification on Ecoli dataset
# ======================================================

library(igraph)
library(ggplot2)
library(reshape2)
library(caret)

# ------------------------------------------------------
# Load TopKGraphs utilities
# ------------------------------------------------------
source("/home/bpfeif/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/calc_BINARY.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bpfeif/GitHub/TopKGraphs/simulations/call_kNN.R")
source("/home/bpfeif/GitHub/TopKGraphs/simulations/call_node2vec.R")

library(mlbench)

data(BreastCancer)

bc <- BreastCancer

# Remove ID column
bc$Id <- NULL

# Remove rows with missing values
bc <- na.omit(bc)

# Features (convert factors to numeric)
X <- as.matrix(
  data.frame(lapply(bc[, -ncol(bc)], function(x) as.numeric(as.character(x))))
)

# Class labels: benign / malignant
labels_true <- as.numeric(as.factor(bc$Class))

n_nodes <- nrow(X)

# Scale features
X <- scale(X)

# ------------------------------------------------------
# kNN graph construction
# ------------------------------------------------------
create_knn_graph <- function(X, k = 5) {
  dist_mat <- as.matrix(dist(X))
  g <- make_empty_graph(n = nrow(X), directed = FALSE)

  for (i in 1:nrow(X)) {
    nn <- order(dist_mat[i, ])[2:(k + 1)]
    for (j in nn) {
      g <- add_edges(g, c(i, j))
    }
  }

  simplify(g)
}

# ------------------------------------------------------
# Safe accuracy computation (FIX)
# ------------------------------------------------------
safe_accuracy <- function(pred, truth) {
  all_levels <- sort(unique(c(pred, truth)))
  cm <- confusionMatrix(
    factor(pred, levels = all_levels),
    factor(truth, levels = all_levels)
  )
  cm$overall["Accuracy"]
}

# ------------------------------------------------------
# Experiment parameters
# ------------------------------------------------------
ks <- c(3, 5, 7, 10)    # kNN graph sizes
setK <- 10             # kNN classifier
n_runs <- 10           # robustness runs

methods <- c(
  "TopKGraphs", "Jaccard", "Dice", "PageRank", "Node2Vec"
)

ACC_array <- array(
  NA,
  dim = c(length(ks), length(methods), n_runs),
  dimnames = list(paste0("k=", ks), methods, NULL)
)

# ------------------------------------------------------
# Main evaluation loop
# ------------------------------------------------------
for (idx in seq_along(ks)) {

  k_graph <- ks[idx]
  cat("Processing kNN graph with k =", k_graph, "\n")

  for (run in 1:n_runs) {

    # ---- robustness noise ----
    X_perturbed <- X + matrix(
      rnorm(n_nodes * ncol(X), sd = 1e-5),
      n_nodes, ncol(X)
    )

    # ---- kNN graph ----
    g <- create_knn_graph(X_perturbed, k = k_graph)

    # ---- TopKGraphs ----
    res_tkg <- topkgraphs(list(g), walk_depth = 30, n_iter = 50)
    dist_tkg <- res_tkg$DIST

    # ---- Jaccard / Dice ----
    jaccard_sim <- similarity(g, method = "jaccard", mode = "all")
    dice_sim <- similarity(g, method = "dice", mode = "all")

    # ---- Personalized PageRank ----
    n <- vcount(g)
    ppr_mat <- sapply(seq_len(n), function(i) {
      pers <- rep(0, n)
      pers[i] <- 1
      page_rank(g, personalized = pers, damping = 0.7)$vector
    })

    # ---- Node2Vec ----
    node2vec_emb <- call_node2vec(g, walk_length = 30, num_walks = 50)
    node_order <- round(as.numeric(node2vec_emb[, 1]))
    ids <- match(1:n_nodes, node_order)
    node2vec_emb <- as.matrix(node2vec_emb[ids, -1])

    # ---- kNN classification ----
    knn_tkg  <- call_kNN_dist(dist_tkg, labels_true, k = setK)
    knn_jac  <- call_kNN_dist(1 - jaccard_sim, labels_true, k = setK)
    knn_dice <- call_kNN_dist(1 - dice_sim, labels_true, k = setK)
    knn_ppr  <- call_kNN_dist(1 - ppr_mat, labels_true, k = setK)
    knn_n2v  <- call_kNN_dist(
      as.matrix(dist(scale(node2vec_emb))),
      labels_true, k = setK
    )

    # ---- Accuracy (SAFE) ----
    ACC_array[idx,"TopKGraphs",run] <- safe_accuracy(knn_tkg[[1]], knn_tkg[[2]])
    ACC_array[idx,"Jaccard",run]   <- safe_accuracy(knn_jac[[1]], knn_jac[[2]])
    ACC_array[idx,"Dice",run]      <- safe_accuracy(knn_dice[[1]], knn_dice[[2]])
    ACC_array[idx,"PageRank",run]  <- safe_accuracy(knn_ppr[[1]], knn_ppr[[2]])
    ACC_array[idx,"Node2Vec",run]  <- safe_accuracy(knn_n2v[[1]], knn_n2v[[2]])
  
  print(ACC_array)
  }
  
}

# ------------------------------------------------------
# Aggregate results
# ------------------------------------------------------
ACC_mean <- apply(ACC_array, c(1, 2), mean)
ACC_sd   <- apply(ACC_array, c(1, 2), sd)

ACC_df <- data.frame(
  k = rep(ks, each = length(methods)),
  Method = rep(methods, times = length(ks)),
  Accuracy = as.vector(t(ACC_mean)),
  SD = as.vector(t(ACC_sd))
)

# ------------------------------------------------------
# Plot
# ------------------------------------------------------
ggplot(ACC_df, aes(x = factor(k), y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymin = Accuracy - SD, ymax = Accuracy + SD),
    width = 0.2, position = position_dodge(0.9)
  ) +
  ylim(0, 1) +
  labs(
    title = "kNN Classification Accuracy (Ecoli)",
    x = "k (kNN graph)",
    y = "Accuracy"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")