# ======================================================
# 0. Libraries
# ======================================================
library(curatedTCGAData)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(igraph)
library(matrixStats)
library(TopKLists)       # TopKGraphs dependencies
library(aricode)         # ARI, NMI

# ======================================================
# 1. Load TCGA BRCA multi-omics
# ======================================================
brca <- curatedTCGAData(
  diseaseCode = "BRCA",
  assays = c("RNASeq2GeneNorm", "Methylation", "miRNASeqGene"),
  dry.run = FALSE,
  version = "2.1.1"
)

# Extract aligned samples for each omics
X_rna  <- t(assay(brca[["BRCA_RNASeq2GeneNorm-20160128"]]))
#X_met  <- t(assay(brca[["BRCA_Methylation"]]))
#X_mir  <- t(assay(brca[["BRCA_miRNASeqGene"]]))

# Find common patients
#common_ids <- Reduce(intersect, list(rownames(X_rna), rownames(X_met), rownames(X_mir)))
#X_rna <- X_rna[common_ids, ]
#X_met <- X_met[common_ids, ]
#X_mir <- X_mir[common_ids, ]

# ======================================================
# 2. Extract PAM50 labels (ground truth)
# ======================================================
clinical <- colData(brca)

# Align X_rna with clinical
X_rna_names = rownames(X_rna)
X_rna_names_clean <- sub("-[0-9][0-9]$", "", X_rna_names)
clinical_names = clinical[,1]

ids = match(X_rna_names_clean, clinical_names)

clinical_clean = clinical[ids,]

labels <- clinical_clean$PAM50.mRNA

valid = which(!is.na(labels))

# Subset omics to valid labels
X_rna_2  = X_rna[valid, ]
labels_2 = labels[valid]

# Remove Normal-like and HER2-enriched
iii = which(labels_2=="Normal-like" | labels_2=="HER2-enriched")
X_rna_2 = X_rna_2[-iii,]
labels_2 = labels_2[-iii]

k <- length(unique(labels_2))
cat("Number of classes:", k, "\n")
table(labels_2)

# ======================================================
# 3. Feature filtering (top variable features)
# ======================================================
top_var <- function(X, p = 1000){
  vars <- rowVars(t(X))
  X[, order(vars, decreasing = TRUE)[1:p]]
}

X_rna_3 <- top_var(X_rna_2, 1000)

# ======================================================
# 4. Build kNN graphs per omics
# ======================================================
build_knn_graph <- function(X, k = 10){
  D <- as.matrix(dist(scale(X)))
  A <- matrix(0, nrow(X), nrow(X))
  for(i in 1:nrow(X)){
    nn <- order(D[i, ])[2:(k+1)]
    A[i, nn] <- 1
  }
  A <- (A + t(A)) > 0
  graph_from_adjacency_matrix(A, mode = "undirected")
}

g_rna <- build_knn_graph(X_rna_3, 20)
graphs <- list(g_rna)

# ======================================================
# 5. Load TopKGraphs functions
# ======================================================
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs_walk.R")
Rcpp::sourceCpp("/home/bastian/GitHub/TopKGraphs/src/walk_with_jaccard_degree_safe.cpp")

# ======================================================
# 6. Run TopKGraphs
# ======================================================
res <- topkgraphs(
  views = graphs,
  walk_depth = 300,
  n_iter = 100,
  do.BORDA = TRUE,
  n_cores=5
)

# ======================================================
# 7. Clustering
# ======================================================
hc <- hclust(as.dist(res$DIST), method = "ward.D2")
cl_topk <- cutree(hc, k = k)

# ======================================================
# 8. Evaluation
# ======================================================
ARI <- ARI(cl_topk, labels_2)
NMI <- NMI(cl_topk, labels_2)
cat("Adjusted Rand Index (ARI):", ARI, "\n")
cat("Normalized Mutual Information (NMI):", NMI, "\n")


P = cbind(cmdscale(res$DIST_raw,2),labels_2)

P_df <- data.frame(
  x = as.numeric(P[, 1]),
  y = as.numeric(P[, 2]),
  label = factor(P[, 3])
)

colnames(P_df) = c("x","y", "label")

library(ggplot2)

ggplot(P_df, aes(x = x, y = y, color = label)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Dimension 1",
    y = "Dimension 2",
    color = "Subtype",
    title = "Patient embedding colored by subtype"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# ======================================================
# 9. Optional: visualize distance matrix
# ======================================================
library(pheatmap)
pheatmap(res$DIST,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "TopKGraphs Patient Similarity",
         show_rownames = FALSE,
         show_colnames = FALSE)
