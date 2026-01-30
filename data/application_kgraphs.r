# ======================================================
# TopKGraphs robustness evaluation on Iris dataset
# ======================================================

library(igraph)
library(ggplot2)
library(reshape2)
library(igraph)
library(ggplot2)
library(reshape2)
library(aricode)  # for ARI/NMI
# Make sure your TopKGraphs functions are loaded: topkgraphs(), etc.


source("/home/bpfeif/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/calc_BINARY.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bpfeif/GitHub/TopKGraphs/simulations/call_kNN.R")
source("/home/bpfeif/GitHub/TopKGraphs/simulations/call_node2vec.R")
# ======================================================
# TopKGraphs + other methods evaluation on Iris dataset
# ======================================================
# ======================================================
# TopKGraphs robustness evaluation on Iris dataset
# ======================================================

library(igraph)
library(ggplot2)
library(reshape2)
library(aricode)  # for ARI/NMI

# Load TopKGraphs functions
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bpfeif/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bpfeif/GitHub/TopKGraphs/simulations/call_node2vec.R")

# -------------------------------
# Load Ecoli dataset from UCI
# -------------------------------
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/ecoli/ecoli.data"

# Column names from the dataset description
cols <- c("SequenceName", "mcg", "gvh", "lip", "chg", "aac", "alm1", "alm2", "Class")

ecoli <- read.table(url, sep="", col.names = cols, stringsAsFactors = FALSE)

# Remove the non-numeric identifier column
X <- as.matrix(ecoli[, 2:8])

# Extract class labels
labels_true <- as.numeric(as.factor(ecoli$Class))

n_nodes <- nrow(X)

# Optional: scale features
X <- scale(X)


# -------------------------------
# 2. Function to create kNN graph
# -------------------------------
create_knn_graph <- function(X, k = 5, metric = "euclidean") {
  dist_mat <- as.matrix(dist(X, method = metric))
  g <- make_empty_graph(n = nrow(X), directed = FALSE)
  for(i in 1:nrow(X)){
    neighbors <- order(dist_mat[i,])[2:(k+1)]
    for(j in neighbors){
      g <- add_edges(g, c(i, j))
    }
  }
  g <- simplify(g)
  return(g)
}

# -------------------------------
# 3. Evaluate multiple kNN graphs with stability
# -------------------------------
ks <- c(3,5,7,10)
methods <- c("TopKGraphs", "Jaccard", "Dice", "PageRank", "Node2Vec")
n_runs <- 10  # number of runs to assess stability

# Store ARI for each run
ARI_array <- array(NA, dim = c(length(ks), length(methods), n_runs),
                   dimnames = list(paste0("k=", ks), methods, NULL))

#set.seed(123)  # reproducibility

for(idx in seq_along(ks)){
  k <- ks[idx]
  cat("Processing k =", k, "\n")
  
  for(run in 1:n_runs){
    # Optional: add small perturbation/noise to X for robustness check
    X_perturbed <- X + matrix(rnorm(n_nodes*ncol(X), mean=0, sd=1e-5), n_nodes, ncol(X))
    
    # ----- kNN graph -----
    g <- create_knn_graph(X_perturbed, k = k)
    
    # ----- TopKGraphs -----
    res_tkg <- topkgraphs(list(g), walk_depth = 30, n_iter = 50,
                          do.BORDA = TRUE, do.TopKSignal = FALSE, do.RRA = FALSE)
    hc_tkg <- hclust(as.dist(res_tkg$DIST), method = "ward.D2")
    cl_tkg <- cutree(hc_tkg, k = length(unique(labels_true)))
    
    # ----- Jaccard & Dice -----
    jaccard_sim <- similarity(g, method = "jaccard", mode = "all")
    dice_sim <- similarity(g, method = "dice", mode = "all")
    
    hc_jaccard <- hclust(as.dist(1-jaccard_sim), method = "ward.D2")
    cl_jaccard <- cutree(hc_jaccard, k = length(unique(labels_true)))
    
    hc_dice <- hclust(as.dist(1-dice_sim), method = "ward.D2")
    cl_dice <- cutree(hc_dice, k = length(unique(labels_true)))
    
    # ----- Personalized PageRank -----
    n <- vcount(g)
    ppr_mat <- sapply(1:n, function(i){
      pers <- rep(0, n)
      pers[i] <- 1
      page_rank(g, personalized = pers, damping=0.7)$vector
    })
    hc_ppr <- hclust(as.dist(1-ppr_mat), method="ward.D2")
    cl_ppr <- cutree(hc_ppr, k = length(unique(labels_true)))
    
    # ----- Node2Vec embeddings -----
    node2vec_emb <- call_node2vec(g, walk_length = 30, num_walks = 50)
    node_order <- round(as.numeric(node2vec_emb[,1]))
    node2vec_emb <- as.matrix(node2vec_emb[,-1])
    hc_n2v <- hclust(dist(scale(node2vec_emb)), method = "ward.D2")
    cl_n2v <- cutree(hc_n2v, k = length(unique(labels_true)))
    ids <- match(1:n_nodes, node_order)
    cl_n2v <- cl_n2v[ids]
    
    # ----- Store ARI -----
    ARI_array[idx,"TopKGraphs",run] <- ARI(labels_true, cl_tkg)
    ARI_array[idx,"Jaccard",run]   <- ARI(labels_true, cl_jaccard)
    ARI_array[idx,"Dice",run]      <- ARI(labels_true, cl_dice)
    ARI_array[idx,"PageRank",run]  <- ARI(labels_true, cl_ppr)
    ARI_array[idx,"Node2Vec",run]  <- ARI(labels_true, cl_n2v)
  
  print(ARI_array)
  }
  
}

# -------------------------------
# 4. Compute mean + SD across runs
# -------------------------------
ARI_mean <- apply(ARI_array, c(1,2), mean)
ARI_sd   <- apply(ARI_array, c(1,2), sd)

ARI_df <- data.frame(
  k = rep(ks, each = length(methods)),
  Method = rep(methods, times = length(ks)),
  ARI_mean = as.vector(t(ARI_mean)),
  ARI_sd   = as.vector(t(ARI_sd))
)
# -------------------------------
# 5. Visualize mean ± SD
# -------------------------------
ggplot(ARI_df, aes(x = factor(k), y = ARI_mean, fill = Method)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin = ARI_mean - ARI_sd, ymax = ARI_mean + ARI_sd),
                width = 0.2, position=position_dodge(0.9)) +
  ylim(0,1) +
  labs(title = "Clustering performance (ARI) across kNN graphs with stability",
       x = "k (kNN graph)", y = "Adjusted Rand Index") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")