# TopKGraphs

# TopKGraphs

**TopKGraphs** computes node affinities in graphs using multiple biased random walks, aggregates partial rankings via Borda counts, and produces a distance matrix suitable for clustering, visualization, or downstream analysis.

---

## Installation

```r
# Install from CRAN (if published)
# install.packages("TopKGraphs")

# Or install latest version from GitHub
# devtools::install_github("pievos101/TopKGraphs")

library(TopKGraphs)
```

## Minimal Example

```r

library(igraph)
library(TopKGraphs)
library(aricode)
library(netUtils)   # For LFR benchmark graphs

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
```


## Run TopKGraphs
```r
res <- topkgraphs(
  views = list(g),
  walk_depth = 50,
  n_iter = 50
)
```

## Baseline: Jaccard similarity
```r
jaccard_sim <- similarity(g, method = "jaccard", mode = "all")
```

## Clustering & Evaluation

```r
k <- length(unique(V(g)$community))

cl_topk <- cutree(hclust(as.dist(res$DIST), method = "ward.D2"), k)
cl_jaccard <- cutree(hclust(as.dist(1 - jaccard_sim), method = "ward.D2"), k)

ari_topk <- ARI(V(g)$community, cl_topk)
ari_jaccard <- ARI(V(g)$community, cl_jaccard)

cat("Adjusted Rand Index (ARI)\n")
cat("-------------------------\n")
cat("TopKGraphs:", round(ari_topk, 3), "\n")
cat("Jaccard:   ", round(ari_jaccard, 3), "\n")
```

## Visualizing Node Affinities
You can visualize the graph embeddings derived from TopKGraphs:

```r
coords <- cmdscale(as.dist(res$DIST_raw), k = 2)
plot(coords, col = V(g)$community, pch = 19,
     main = "2D Embedding from TopKGraphs Distances")
```