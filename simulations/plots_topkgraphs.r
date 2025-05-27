# Convert to graph object
graph_fused <- graph_from_adjacency_matrix(
  res,
  mode = "directed",
  weighted = TRUE,
  diag = FALSE
)

# Plot the graph with edge weights as labels
plot(
  graph_fused,
  edge.label = E(graph_fused)$weight,  # This shows weights on the edges
  edge.arrow.size = 0.5,               # Make arrows smaller
  vertex.label.cex = 1,              # Vertex label size
  edge.label.cex = 1,                  # Edge label size
  edge.label.color = "darkred"         # Edge label color
)