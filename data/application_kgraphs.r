# ======================================================
# 0. Libraries
# ======================================================
library(igraph)
library(ggplot2)
library(reshape2)
library(aricode)   # ARI, NMI, AMI

# ======================================================
# 1. Load TopKGraphs functions
# ======================================================
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bastian/GitHub/TopKGraphs/simulations/call_node2vec.R")

# -------------------------------
# Load Breast Cancer dataset (UCI)
# -------------------------------
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"

cols <- c(
  "ID", "Diagnosis",
  paste0("Feature", 1:30)
)

bc <- read.table(url, sep = ",", col.names = cols, stringsAsFactors = FALSE)

# Features (exclude ID and Diagnosis)
X <- as.matrix(bc[, -c(1, 2)])

# Ground-truth labels (M = malignant, B = benign)
labels_true <- as.numeric(as.factor(bc$Diagnosis))

# Number of samples
n_nodes <- nrow(X)

# Scale features
X <- scale(X)

# Number of clusters
n_clusters <- length(unique(labels_true))


# ======================================================
# 3. kNN graph constructor
# ======================================================
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

# ======================================================
# 4. External metrics helper
# ======================================================
compute_external_metrics <- function(labels_true, labels_pred) {
  c(
    ARI = ARI(labels_true, labels_pred),
    NMI = NMI(labels_true, labels_pred),
    AMI = AMI(labels_true, labels_pred)
  )
}

# ======================================================
# 5. Experiment setup
# ======================================================
ks <- c(5, 7, 10)
methods <- c("TopKGraphs", "Jaccard", "Dice", "PageRank", "Node2Vec")
metrics <- c("ARI", "NMI", "AMI")
n_runs <- 10

METRIC_array <- array(
  NA,
  dim = c(length(ks), length(methods), n_runs, length(metrics)),
  dimnames = list(
    paste0("k=", ks),
    methods,
    NULL,
    metrics
  )
)

# ======================================================
# 6. Main loop
# ======================================================
for (idx in seq_along(ks)) {

  k <- ks[idx]
  cat("Processing k =", k, "\n")

  for (run in 1:n_runs) {

    # small perturbation for robustness
    X_perturbed <- X + matrix(
      rnorm(n_nodes * ncol(X), sd = 1e-5),
      n_nodes, ncol(X)
    )

    # ---- kNN graph ----
    g <- create_knn_graph(X_perturbed, k = k)

    # ---- TopKGraphs ----
    res_tkg <- topkgraphs(
      list(g),
      walk_depth = 30,
      n_iter = 50,
      do.BORDA = TRUE,
      do.TopKSignal = FALSE,
      do.RRA = FALSE
    )

    hc_tkg <- hclust(as.dist(res_tkg$DIST), method = "ward.D2")
    cl_tkg <- cutree(hc_tkg, k = n_clusters)

    METRIC_array[idx, "TopKGraphs", run, ] <-
      compute_external_metrics(labels_true, cl_tkg)

    # ---- Jaccard / Dice ----
    jaccard_sim <- similarity(g, method = "jaccard", mode = "all")
    dice_sim <- similarity(g, method = "dice", mode = "all")

    cl_jaccard <- cutree(
      hclust(as.dist(1 - jaccard_sim), method = "ward.D2"),
      k = n_clusters
    )
    cl_dice <- cutree(
      hclust(as.dist(1 - dice_sim), method = "ward.D2"),
      k = n_clusters
    )

    METRIC_array[idx, "Jaccard", run, ] <-
      compute_external_metrics(labels_true, cl_jaccard)
    METRIC_array[idx, "Dice", run, ] <-
      compute_external_metrics(labels_true, cl_dice)

    # ---- Personalized PageRank ----
    n <- vcount(g)
    ppr_mat <- sapply(seq_len(n), function(i) {
      pers <- rep(0, n)
      pers[i] <- 1
      page_rank(g, personalized = pers, damping = 0.7)$vector
    })

    cl_ppr <- cutree(
      hclust(as.dist(1 - ppr_mat), method = "ward.D2"),
      k = n_clusters
    )

    METRIC_array[idx, "PageRank", run, ] <-
      compute_external_metrics(labels_true, cl_ppr)

    # ---- Node2Vec ----
    node2vec_emb <- call_node2vec(g, walk_length = 30, num_walks = 50)
    node_order <- round(as.numeric(node2vec_emb[, 1]))
    emb <- as.matrix(node2vec_emb[, -1])

    hc_n2v <- hclust(dist(scale(emb)), method = "ward.D2")
    cl_n2v <- cutree(hc_n2v, k = n_clusters)

    ids <- match(1:n_nodes, node_order)
    cl_n2v <- cl_n2v[ids]

    METRIC_array[idx, "Node2Vec", run, ] <-
      compute_external_metrics(labels_true, cl_n2v)

      print(METRIC_array)
  }
}

# ======================================================
# 7. Aggregate results
# ======================================================
METRIC_mean <- apply(METRIC_array, c(1, 2, 4), mean, na.rm = TRUE)
METRIC_sd   <- apply(METRIC_array, c(1, 2, 4), sd, na.rm = TRUE)

METRIC_df <- melt(METRIC_mean)
colnames(METRIC_df) <- c("k", "Method", "Metric", "Mean")
METRIC_df$SD <- as.vector(METRIC_sd)

# ======================================================
# 8. Plot
# ======================================================
ggplot(METRIC_df, aes(x = k, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Clustering performance on BreastCancer dataset",
    x = "k (kNN graph)",
    y = "Score"
  )
