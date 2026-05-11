#' Temporal Graphs from longitudinal data
#'
#'
#' @param wide Longitudinal data in Wide Format
#' @param k_similarity within time-slice edges
#'
#' @return A list containing the graph, labels, data
#'
#' @export

build_temporal_graph <- function(wide, k_similarity = 10) {

  require(igraph)
  wide <- as.data.frame(wide)

  # -------------------------------------------------
  # STEP 1 — feature matrix
  # -------------------------------------------------
  feature_cols <- grep("^y\\.", names(wide), value = TRUE)

  X <- as.matrix(wide[, feature_cols])

  # -------------------------------------------------
  # STEP 2 — stable node ids
  # -------------------------------------------------
  wide$node_id <- seq_len(nrow(wide)) - 1

  # -------------------------------------------------
  # STEP 3 — patient-level labels
  # -------------------------------------------------
  patient_labels <- unique(wide[, c("subject", "cluster")])

  patient_labels <- patient_labels[
    order(patient_labels$subject),
  ]

  # -------------------------------------------------
  # STEP 4 — create graph + nodes
  # -------------------------------------------------
  g <- make_empty_graph(directed = FALSE)

  g <- add_vertices(
    g,
    nv = nrow(wide),
    name = as.character(wide$node_id),
    subject = wide$subject,
    time = wide$time,
    cluster = wide$cluster
  )

  # feature vectors as vertex attributes
  V(g)$features <- lapply(seq_len(nrow(X)), function(i) {
    as.numeric(X[i, ])
  })

  # -------------------------------------------------
  # STEP 5 — temporal edges
  # -------------------------------------------------
  temporal_edges <- c()
  temporal_types <- c()

  subjects <- unique(wide$subject)

  for (s in subjects) {

    sub <- wide[wide$subject == s, ]
    sub <- sub[order(sub$time), ]

    if (nrow(sub) < 2)
      next

    for (i in 1:(nrow(sub) - 1)) {

      from <- as.character(sub$node_id[i])
      to   <- as.character(sub$node_id[i + 1])

      temporal_edges <- c(
        temporal_edges,
        from,
        to
      )

      temporal_types <- c(
        temporal_types,
        "temporal"
      )
    }
  }

  if (length(temporal_edges) > 0) {

    g <- add_edges(
      g,
      temporal_edges,
      attr = list(
        edge_type = temporal_types
      )
    )
  }

  # -------------------------------------------------
  # STEP 6 — similarity edges within time slices
  # -------------------------------------------------
  similarity_edges <- c()
  similarity_types <- c()
  similarity_times <- c()

  unique_times <- sort(unique(wide$time))

  for (t in unique_times) {

    slice_idx <- which(wide$time == t)

    if (length(slice_idx) <= k_similarity)
      next

    X_slice <- X[slice_idx, , drop = FALSE]

    # pairwise Euclidean distances
    dist_mat <- as.matrix(dist(X_slice))

    for (i_local in seq_along(slice_idx)) {

      i_global <- slice_idx[i_local]

      # nearest neighbors excluding self
      neighbors <- order(dist_mat[i_local, ])[2:(k_similarity + 1)]

      for (j_local in neighbors) {

        j_global <- slice_idx[j_local]

        from <- as.character(wide$node_id[i_global])
        to   <- as.character(wide$node_id[j_global])

        similarity_edges <- c(
          similarity_edges,
          from,
          to
        )

        similarity_types <- c(
          similarity_types,
          "similarity"
        )

        similarity_times <- c(
          similarity_times,
          t
        )
      }
    }
  }

  if (length(similarity_edges) > 0) {

    g <- add_edges(
      g,
      similarity_edges,
      attr = list(
        edge_type = similarity_types,
        time = similarity_times
      )
    )
  }

## INFO
#V(g)$name
#V(g)$subject
#V(g)$time

  return(list(
    graph = g,
    patient_labels = patient_labels,
    wide_data = wide
  ))
}