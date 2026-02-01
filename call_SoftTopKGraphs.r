library(reticulate)

# -----------------------------
# 1️⃣ Create or use a dedicated virtual environment
# -----------------------------
env_name <- "topkgraphs_env"

# Check if the env exists; if not, create it
if (!virtualenv_exists(env_name)) {
  cat("Creating virtual environment...\n")
  virtualenv_create(env_name, python = "/usr/bin/python3")
}

# Activate the virtual environment
use_virtualenv(env_name, required = TRUE)

# -----------------------------
# 2️⃣ Install PyTorch if missing
# -----------------------------
if (!py_module_available("torch")) {
  cat("Installing PyTorch...\n")
  py_install("torch", envname = env_name, pip = TRUE)
}

# -----------------------------
# 3️⃣ Import your Python module
# -----------------------------
soft_topkgraphs <- import_from_path(
  "SoftTopKGraphs", 
  path = "/home/bpfeif/GitHub/TopKGraphs/"
)

# -----------------------------
# 4️⃣ Quick test
# -----------------------------
# Example: create a small graph in Python via networkx
nx <- import("networkx")
G <- nx$cycle_graph(5L)  # note the L to force integer in Python

# Call the SoftTopKGraphs class
model <- soft_topkgraphs$SoftTopKGraphs(G, walk_length = 5L)

# Generate embeddings
embeddings <- model()$detach()$numpy()
print(embeddings)