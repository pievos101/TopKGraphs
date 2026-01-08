#
# python3 -m venv pyg-env
# source pyg-env/bin/activate
# pip install "numpy<2.0"
# 
# deactivate
# ipython3

#######################
# Citeseer
#######################
import numpy as np
from torch_geometric.datasets import Planetoid

dataset = Planetoid(root="data/Citeseer", name="Citeseer")

data = dataset[0]

# Get node labels (shape: [num_nodes])
labels = data.y

# Convert to NumPy
labels_np = labels.cpu().numpy()

# OPTIONAL: convert to 1-based labels for R
labels_np = labels_np + 1

# Save to CSV (one label per line)
np.savetxt(
    "citeseer_node_labels.csv",
    labels_np,
    fmt="%d"
)

edge_index = data.edge_index

edge_matrix = edge_index.t().contiguous()

# Convert to NumPy
edge_matrix_np = edge_matrix.cpu().numpy()

# OPTIONAL: convert to 1-based indexing for R
edge_matrix_np = edge_matrix_np + 1

# Save as CSV
np.savetxt(
    "citeseer_edge_matrix.csv",
    edge_matrix_np,
    fmt="%d",
    delimiter=","
)
#################################################################