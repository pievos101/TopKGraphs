## TopKGraphs: minimal example
## ---------------------------

library(igraph)
library(TopKGraphs)
library(aricode)
library(netUtils)   # <-- For LFR benchmark graphs

# ----------------------------
# Generate a small benchmark graph
# ----------------------------
g <- sample_lfr(
  n = 60,
  average_degree = 6,
  max_degree = 12,
  min_community = 10,
  max_community = 30,
  mu = 0.2,
  tau1 = 2,
  tau2 = 1.2
)

# Ground truth communities
V(g)$community <- V(g)$membership
g <- delete_vertices(g, V(g)[degree(g) == 0])

# ----------------------------
# Run TopKGraphs
# ----------------------------
res <- topkgraphs(
  views = list(g),
  walk_depth = 50,
  n_iter = 50
)

# Convert distance to similarity
topk_sim <- 1 - res$DIST

# ----------------------------
# Baseline: Jaccard similarity
# ----------------------------
jaccard_sim <- similarity(g, method = "jaccard", mode = "all")

# ----------------------------
# Clustering
# ----------------------------
k <- length(unique(V(g)$community))

cl_topk <- cutree(
  hclust(as.dist(res$DIST), method = "ward.D2"),
  k
)

cl_jaccard <- cutree(
  hclust(as.dist(1 - jaccard_sim), method = "ward.D2"),
  k
)

# ----------------------------
# Evaluation
# ----------------------------
ari_topk <- ARI(V(g)$community, cl_topk)
ari_jaccard <- ARI(V(g)$community, cl_jaccard)

cat("Adjusted Rand Index (ARI)\n")
cat("-------------------------\n")
cat("TopKGraphs:", round(ari_topk, 3), "\n")
cat("Jaccard:   ", round(ari_jaccard, 3), "\n")