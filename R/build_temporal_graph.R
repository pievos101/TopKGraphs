#' Temporal Graphs from longitudinal data
#'
#'
#' @param wide Longitudinal data in Wide Format
#' @param k_similarity within time-slice edges
#'
#' @return A list containing the graph, labels, data
#'
#' @export

build_temporal_graph <- function(wide, k_similarity = 10, n_subspaces = 5, subspace_ratio = 0.7) {

  library(igraph)

  wide <- as.data.frame(wide)

  # -------------------------------------------------
  # STEP 1 — feature matrix
  # -------------------------------------------------
  feature_cols <- grep("^y\\.", names(wide), value = TRUE)
  X_raw <- as.matrix(wide[, feature_cols])

  # -------------------------------------------------
  # STEP 1.1 — ROBUST FEATURE NORMALIZATION
  # (critical: removes scale sensitivity + reduces noise impact)
  # -------------------------------------------------
  X <- scale(X_raw)

  # clip extreme values (winsorization)
  clip_val <- 4
  X[X > clip_val] <- clip_val
  X[X < -clip_val] <- -clip_val

  n <- nrow(X)

  # -------------------------------------------------
  # STEP 2 — stable node ids
  # -------------------------------------------------
  wide$node_id <- seq_len(n) - 1

  # -------------------------------------------------
  # STEP 3 — patient-level labels
  # -------------------------------------------------
  patient_labels <- unique(wide[, c("subject", "cluster")])
  patient_labels <- patient_labels[order(patient_labels$subject), ]

  # -------------------------------------------------
  # STEP 4 — create graph
  # -------------------------------------------------
  g <- make_empty_graph(directed = FALSE)

  g <- add_vertices(
    g,
    nv = n,
    name = as.character(wide$node_id),
    subject = wide$subject,
    time = wide$time,
    cluster = wide$cluster
  )

  V(g)$features <- lapply(seq_len(n), function(i) as.numeric(X[i, ]))

  # -------------------------------------------------
  # STEP 5 — temporal edges (UNCHANGED)
  # -------------------------------------------------
  temporal_edges <- c()
  temporal_types <- c()

  subjects <- unique(wide$subject)

  for (s in subjects) {

    sub <- wide[wide$subject == s, ]
    sub <- sub[order(sub$time), ]

    if (nrow(sub) < 2) next

    for (i in 1:(nrow(sub) - 1)) {

      temporal_edges <- c(
        temporal_edges,
        as.character(sub$node_id[i]),
        as.character(sub$node_id[i + 1])
      )

      temporal_types <- c(temporal_types, "temporal")
    }
  }

  if (length(temporal_edges) > 0) {
    g <- add_edges(g, temporal_edges, attr = list(edge_type = temporal_types))
  }

  # -------------------------------------------------
  # STEP 6 — ROBUST kNN SIMILARITY GRAPH
  # -------------------------------------------------
  similarity_edges <- c()
  similarity_types <- c()

  p <- ncol(X)
  subspace_size <- max(2, floor(subspace_ratio * p))

  unique_times <- sort(unique(wide$time))

  for (t in unique_times) {

    idx <- which(wide$time == t)

    if (length(idx) <= k_similarity) next

    X_slice <- X[idx, , drop = FALSE]

    # -------------------------------------------------
    # ENSEMBLE kNN (KEY FIX)
    # -------------------------------------------------
    agg_dist <- matrix(0, nrow = nrow(X_slice), ncol = nrow(X_slice))

    for (b in 1:n_subspaces) {

      feat_idx <- sample(1:p, subspace_size, replace = FALSE)
      X_sub <- X_slice[, feat_idx, drop = FALSE]

      d <- as.matrix(dist(X_sub))
      agg_dist <- agg_dist + rank(d)
    }

    agg_dist <- agg_dist / n_subspaces

    # -------------------------------------------------
    # kNN from aggregated ranks (robust to noise dims)
    # -------------------------------------------------
    for (i in 1:nrow(agg_dist)) {

      neighbors <- order(agg_dist[i, ])[2:(k_similarity + 1)]

      for (j in neighbors) {

        similarity_edges <- c(
          similarity_edges,
          as.character(wide$node_id[idx[i]]),
          as.character(wide$node_id[idx[j]])
        )

        similarity_types <- c(similarity_types, "similarity")
      }
    }
  }

  if (length(similarity_edges) > 0) {
    g <- add_edges(g, similarity_edges, attr = list(edge_type = similarity_types))
  }

  # -------------------------------------------------
  # RETURN
  # -------------------------------------------------
  return(list(
    graph = g,
    patient_labels = patient_labels,
    wide_data = wide
  ))
}