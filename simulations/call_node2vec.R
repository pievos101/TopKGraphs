library(igraph)
library(reticulate)

# Use appropriate Python environment
use_virtualenv("r-reticulate", required = FALSE)  # or use_python("/usr/bin/python3")

# Ensure required Python packages are available
py_install(c("node2vec", "networkx", "pandas"), pip = TRUE)

# TODO -- make walk_length and num_walks params of the function
call_node2vec <- function(g, walk_length=50, num_walks=500){


walk_length = as.integer(walk_length)
num_walks = as.integer(num_walks)

# Create igraph graph
# g <- make_ring(10) %>% add_edges(c(1,5, 2,6, 3,7))

# Convert to edge list
edge_list <- as.data.frame(get.edgelist(g))
colnames(edge_list) <- c("source", "target")

# Load pandas and pass as raw Python object
pd <- import("pandas")
edges_dict <- r_to_py(edge_list, convert = FALSE)
py$edges_dict <- edges_dict
py$walk_length <- walk_length
py$num_walks <- num_walks

# Use pandas.DataFrame in Python explicitly
py_run_string("
import pandas as pd
from node2vec import Node2Vec
import networkx as nx

edges = pd.DataFrame(edges_dict)
edge_list = edges.to_numpy().tolist()

G = nx.Graph()
G.add_edges_from(edge_list)

node2vec = Node2Vec(G, dimensions=128, walk_length=walk_length, 
        num_walks=num_walks, workers=1)
model = node2vec.fit(window=10, min_count=1, batch_words=4)
embeddings = {str(node): model.wv[str(node)] for node in G.nodes()}
")

# Back to R
embedding_list <- py$embeddings
embedding_matrix <- do.call(rbind, embedding_list)
embedding_df <- data.frame(node = names(embedding_list), embedding_matrix)

return(embedding_df)

}