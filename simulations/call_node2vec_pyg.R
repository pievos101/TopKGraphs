library(igraph)
library(reticulate)

use_virtualenv("r-reticulate", required = FALSE)

# Ensure Python deps
py_install(
  c("torch", "torch-geometric", "torch-cluster", "torch-scatter", "torch-sparse"),
  pip = TRUE
)

call_node2vec_pyg <- function(g,
                             walk_length = 20,
                             num_walks = 10,
                             dimensions = 128,
                             window_size = 10,
                             epochs = 100,
                             lr = 0.01) {

  walk_length <- as.integer(walk_length)
  num_walks <- as.integer(num_walks)
  dimensions <- as.integer(dimensions)
  window_size <- as.integer(window_size)
  epochs <- as.integer(epochs)

  # ---- Convert igraph to edge list ----
  edge_list <- as.data.frame(get.edgelist(g))
  colnames(edge_list) <- c("source", "target")

  # Map node names to 0..N-1 (required by PyG)
  nodes <- sort(unique(c(edge_list$source, edge_list$target)))
  node_map <- setNames(seq_along(nodes) - 1, nodes)

  edge_list$source <- node_map[edge_list$source]
  edge_list$target <- node_map[edge_list$target]

  py$edges <- r_to_py(edge_list)
  py$num_nodes <- length(nodes)
  py$walk_length <- walk_length
  py$num_walks <- num_walks
  py$dimensions <- dimensions
  py$window_size <- window_size
  py$epochs <- epochs
  py$lr <- lr

  # ---- Python: PyTorch Geometric Node2Vec ----
  py_run_string("
import torch
from torch_geometric.nn import Node2Vec

# Edge index tensor [2, num_edges]
edge_index = torch.tensor(
    edges[['source', 'target']].values.T,
    dtype=torch.long
)

device = 'cuda' if torch.cuda.is_available() else 'cpu'

model = Node2Vec(
    edge_index,
    embedding_dim=dimensions,
    walk_length=walk_length,
    context_size=window_size,
    walks_per_node=num_walks,
    num_nodes=num_nodes,
    sparse=True
).to(device)

optimizer = torch.optim.SparseAdam(model.parameters(), lr=lr)

def train():
    model.train()
    total_loss = 0
    for pos_rw, neg_rw in model.loader(batch_size=128, shuffle=True):
        optimizer.zero_grad()
        loss = model.loss(pos_rw.to(device), neg_rw.to(device))
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    return total_loss

for epoch in range(epochs):
    train()

embeddings = model.embedding.weight.detach().cpu().numpy()
")

  # ---- Back to R ----
  embedding_matrix <- py$embeddings
  rownames(embedding_matrix) <- nodes

  #return(embedding_matrix)

  embedding_df <- data.frame(
    node = nodes,
    embedding_matrix,
    row.names = NULL
  )

  return(embedding_df)
}