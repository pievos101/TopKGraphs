# TopKSignal application
# Compare also to https://github.com/berenslab/graph-ne-paper
# Simulate synthetic graphs
# TopKGraphs
library(TopKLists)
library(igraph)
library(fastcluster) 
library(cluster) 
#library(node2vec)


source("/home/bastian/GitHub/TopKGraphs/simulations/sim.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_BINARY.R")
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs.R")
source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs_walk.R")
source("/home/bastian/GitHub/TopKGraphs/simulations/call_node2vec.R")
source("/home/bastian/GitHub/TopKGraphs/simulations/call_kNN.R")


# Sizes of the communities
sizes <- c(10, 10, 10)
n_nodes = sum(sizes)

# Connection probability matrix (3x3)
intra = 0.50 # 0.50 is baseline
inter = 0.05 # 0.05 is baseline

pref.matrix <- matrix(c(
    intra, inter, inter,
    inter, intra, inter,
    inter, inter, intra
), nrow = 3, byrow = TRUE)

# Generate the SBM graph
g <- sample_sbm(sum(sizes), pref.matrix, block.sizes = sizes)

# Assign community membership
V(g)$community <- rep(1:length(sizes), times = sizes)

# Remove isolated nodes
g <- delete_vertices(g, V(g)[degree(g) == 0])
n_nodes = length(V(g))


# Get the rank Matrix for a node

res = topkgraphs_walk(list(g), start_node=20, walk_depth=20, n_iter=30)

RM = res$RM 
colnames(RM) = paste("r",1:ncol(RM), sep="")
rownames(RM) = paste("v",1:nrow(RM), sep="")


IN = matrix(NaN, nrow(RM), ncol(RM))
colnames(IN) = paste("r",1:ncol(RM), sep="")

node_names = sort(unique(as.vector(RM)))


IN = apply(RM, 2, function(x){match(node_names,x)})
rownames(IN) = node_names

library(TopKSignal)
library(gurobi)

estimatedSignal <- estimateTheta(R.input = IN, 
         num.boot = 50, b = 0.1, solver = "gurobi", 
         type = "restrictedQuadratic", bootstrap.type = "classic.bootstrap",
         nCore = 1)

estimatedSignal$estimation

