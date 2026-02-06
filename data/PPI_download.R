# -------------------------------
# 2. Download STRING files
# -------------------------------
# PPI edges
edges_url <- "https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz"
edges_file <- "9606.protein.links.v11.5.txt.gz"
if(!file.exists(edges_file)) download.file(edges_url, edges_file)
edges <- read.table(edges_file, header=TRUE, sep=" ")

# Protein info (maps ENSP -> gene name)
info_url <- "https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz"
info_file <- "9606.protein.info.v11.5.txt.gz"
if(!file.exists(info_file)) download.file(info_url, info_file)
info <- read.table(info_file, sep="\t")

colnames(info) = c("protein_id","gene","3","4")

# Mapping 
edges <- merge(
  edges,
  info[, c("protein_id", "gene")],
  by.x = "protein1",
  by.y = "protein_id",
  all.x = TRUE
)

colnames(edges)[ncol(edges)] <- "gene1"

edges <- merge(
  edges,
  info[, c("protein_id", "gene")],
  by.x = "protein2",
  by.y = "protein_id",
  all.x = TRUE
)

colnames(edges)[ncol(edges)] <- "gene2"

# Remove NaNs
edges <- edges[!is.na(edges$gene1) & !is.na(edges$gene2), ]


PPI = edges[,c("gene1","gene2","combined_score")]

# filter 
ids = which(PPI[,"combined_score"]>=950)
PPI_sub = PPI[ids,]

##################################
# Get the diseases
##################################
# Zenodo DisGeNET gene-disease associations

url <- "https://zenodo.org/records/48426/files/disgenet-v1.0.zip"
zipfile <- "disgenet-v1.0.zip"

if(!file.exists(zipfile)){
    download.file(url, zipfile, mode="wb")
}

unzip(zipfile, exdir="disgenet_v1.0")

# Inspect files
files <- list.files("disgenet_v1.0", full.names=TRUE)
print(files)

# Load gene-disease associations
disease <- read.table(
  "disgenet_v1.0/dhimmel-disgenet-a636cc0/data/consolidated.tsv",
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = TRUE,
  comment.char = "")

# Create the IGRAPH
library(igraph)
g <- graph_from_data_frame(PPI_sub[,c("gene1","gene2")], directed=FALSE)
adj_matrix <- as.matrix(as_adjacency_matrix(g, sparse=FALSE))
cat("PPI network:", vcount(g), "nodes,", ecount(g), "edges\n")

# Connected components
comps <- components(g)

# ID of the largest component
largest_comp_id <- which.max(comps$csize)

# Vertices belonging to the largest component
largest_vertices <- V(g)[comps$membership == largest_comp_id]

# Induce subgraph
g_largest <- induced_subgraph(g, largest_vertices)

cat("Largest component:",
    vcount(g_largest), "nodes,",
    ecount(g_largest), "edges\n")

# Now set the disease communities
ids = match(V(g_largest)$name, disease[,"geneSymbol"])

V(g_largest)$community = disease[ids,"doid_name"]

g_largest = delete_vertices(g_largest, which(is.na(V(g_largest)$community)))

## Again largest component

# Connected components
comps <- components(g_largest)

# ID of the largest component
largest_comp_id <- which.max(comps$csize)

# Vertices belonging to the largest component
largest_vertices <- V(g_largest)[comps$membership == largest_comp_id]

# Induce subgraph
g_largest <- induced_subgraph(g_largest, largest_vertices)

cat("Largest component:",
    vcount(g_largest), "nodes,",
    ecount(g_largest), "edges\n")

# Delete small communities 
ids = which(table(V(g_largest)$community)<50)
nn  = names(table(V(g_largest)$community))[ids]

idx_remove <- which(is.element(V(g_largest)$community, nn))

g_largest = delete_vertices(g_largest, idx_remove)

## Again largest component

# Connected components
comps <- components(g_largest)

# ID of the largest component
largest_comp_id <- which.max(comps$csize)

# Vertices belonging to the largest component
largest_vertices <- V(g_largest)[comps$membership == largest_comp_id]

# Induce subgraph
g_largest <- induced_subgraph(g_largest, largest_vertices)

cat("Largest component:",
    vcount(g_largest), "nodes,",
    ecount(g_largest), "edges\n")


# Now delete ****prostata cancer***** -- too many genes
idx_remove = which(V(g_largest)$community == "prostate cancer")

g_largest = delete_vertices(g_largest, idx_remove)


## Again largest component

# Connected components
comps <- components(g_largest)

# ID of the largest component
largest_comp_id <- which.max(comps$csize)

# Vertices belonging to the largest component
largest_vertices <- V(g_largest)[comps$membership == largest_comp_id]

# Induce subgraph
g_largest <- induced_subgraph(g_largest, largest_vertices)

cat("Largest component:",
    vcount(g_largest), "nodes,",
    ecount(g_largest), "edges\n")


save(g_largest, file="PPI_sub.RD")