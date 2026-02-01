# ======================================================
# TopKGraphs robustness evaluation using multiple metrics
# ======================================================

# -------------------------------
# Load libraries
# -------------------------------
library(igraph)
library(ggplot2)
library(reshape2)
library(cluster)       # silhouette()
library(clusterCrit)   # internal cluster quality metrics

# -------------------------------
# Load TopKGraphs functions
# -------------------------------
source("/home/bpfeif/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/calc_BINARY.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bpfeif/GitHub/TopKGraphs/simulations/call_node2vec.R")

# -------------------------------
# Load Breast Cancer dataset (UCI)
# -------------------------------
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"

# Column names from UCI description
cols <- c(
  "ID",
  "Diagnosis",
  paste0("Feature", 1:30)
)

bc <- read.table(
  url,
  sep = ",",
  col.names = cols,
  stringsAsFactors = FALSE
)

# Features (exclude ID and diagnosis)
X <- as.matrix(bc[, -c(1, 2)])

# Ground-truth class labels (M / B)
labels_true <- as.numeric(as.factor(bc$Diagnosis))

# Number of samples
n_nodes <- nrow(X)

# Scale features
X <- scale(X)

# -------------------------------
# kNN graph constructor
# -------------------------------
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

# -------------------------------
# Metric helpers
# -------------------------------
mean_silhouette <- function(dist_mat, labels) {
  sil <- silhouette(labels, as.dist(dist_mat))
  mean(sil[, "sil_width"])
}

calinski_harabasz <- function(dist_mat, labels) {
  crit <- intCriteria(as.matrix(dist_mat), as.integer(labels), c("Calinski_Harabasz"))
  crit$calinski_harabasz
}

davies_bouldin <- function(dist_mat, labels) {
  crit <- intCriteria(as.matrix(dist_mat), as.integer(labels), c("Davies_Bouldin"))
  crit$davies_bouldin
}

# -------------------------------
# Experiment setup
# -------------------------------
ks <- c(5, 7, 10)
methods <- c("TopKGraphs", "Jaccard", "Dice", "PageRank", "Node2Vec")
n_runs <- 10
metrics <- c("Silhouette", "CalinskiHarabasz", "DaviesBouldin")

# 4D array: k x method x run x metric
RESULT_array <- array(
  NA,
  dim = c(length(ks), length(methods), n_runs, length(metrics)),
  dimnames = list(paste0("k=", ks), methods, NULL, metrics)
)

# -------------------------------
# Main loop
# -------------------------------
for (idx in seq_along(ks)) {

  k <- ks[idx]
  cat("Processing k =", k, "\n")

  for (run in 1:n_runs) {

    # Small perturbation for robustness
    X_perturbed <- X + matrix(rnorm(n_nodes * ncol(X), sd = 1e-5),
                              n_nodes, ncol(X))

    # ---- kNN graph ----
    g <- create_knn_graph(X_perturbed, k = k)

    # ---- TopKGraphs ----
    res_tkg <- topkgraphs(
      list(g),
      walk_depth = 50,
      n_iter = 50,
      do.BORDA = TRUE,
      do.TopKSignal = FALSE,
      do.RRA = FALSE
    )
    dist_tkg = res_tkg$DIST
    
    # ---- Jaccard / Dice ----
    jaccard_sim <- similarity(g, method = "jaccard", mode = "all")
    dice_sim <- similarity(g, method = "dice", mode = "all")
    dist_jaccard <- 1 - jaccard_sim
    dist_dice <- 1 - dice_sim

    # ---- Personalized PageRank ----
    n <- vcount(g)
    ppr_mat <- sapply(seq_len(n), function(i) {
      pers <- rep(0, n)
      pers[i] <- 1
      page_rank(g, personalized = pers, damping = 0.7)$vector
    })
    dist_ppr <- 1 - ppr_mat

    # ---- Node2Vec ----
    node2vec_emb <- call_node2vec(g, walk_length = 50, num_walks = 50)
    node_order <- round(as.numeric(node2vec_emb[, 1]))
    ids <- match(1:n_nodes, node_order)
    node2vec_emb <- as.matrix(node2vec_emb[ids, -1])
    dist_n2v <- as.matrix(dist(scale(node2vec_emb)))

    # ---- Compute metrics ----
    dist_list <- list(
      TopKGraphs = dist_tkg,
      Jaccard = dist_jaccard,
      Dice = dist_dice,
      PageRank = dist_ppr,
      Node2Vec = dist_n2v
    )

    for (m in methods) {
      RESULT_array[idx, m, run, "Silhouette"] <- mean_silhouette(dist_list[[m]], labels_true)
      RESULT_array[idx, m, run, "CalinskiHarabasz"] <- calinski_harabasz(dist_list[[m]], labels_true)
      RESULT_array[idx, m, run, "DaviesBouldin"] <- davies_bouldin(dist_list[[m]], labels_true)
    }

    print(RESULT_array[idx, , run, ])
  }
}

# -------------------------------
# Aggregate results (mean + SD across runs)
# -------------------------------
RESULT_mean <- apply(RESULT_array, c(1,2,4), mean)
RESULT_sd   <- apply(RESULT_array, c(1,2,4), sd)

RESULT_df <- reshape2::melt(RESULT_mean)
colnames(RESULT_df) <- c("k", "Method", "Metric", "Mean")
RESULT_df$SD <- reshape2::melt(RESULT_sd)$value

# -------------------------------
# Plot
# -------------------------------
ggplot(RESULT_df, aes(x = factor(k), y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(
    title = "Cluster quality metrics vs ground truth",
    x = "k (kNN graph)",
    y = "Metric value"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")