# TopKGraphs
library(TopKLists)
library(igraph)
library(fastcluster) 
library(cluster) 

source("/home/bastian/GitHub/TopKGraphs/simulations/sim.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_BINARY.R")

# Varianz 
my_var <- 0.1
VVV <- c(0.1,0.5,1,5,10,15,20)
this_method = "ward.D"
max.k = 7
fix.k = NaN

omics  <- sim1(FALSE, my_var, mode="HC")

omics_binary = list()

for (xx in 1:length(omics)){

	omic  <- omics[[xx]]	

	if(is.na(fix.k)){
		sil   <- calc.SIL(dist(omic), max.k, method=this_method)
		id    <- which.max(sil)
		k     <- as.numeric(names(sil)[id])
	}else{
		k <- fix.k
	}

	hc    <- hclust(dist(omic), method=this_method)
	cl    <- cutree(hc, k)

	mat   <- calc.BINARY(cl)	
        
	omics_binary[[xx]] <- mat 
       
}

source("/home/bastian/GitHub/TopKGraphs/R/topkgraphs.R")

res = topkgraphs(omics_binary)

hc = hclust(as.dist(res$DIST), method="ward.D")


# BASIC PLOT 
plot(hc)



#### ADVANCED PLOTS
###############################################


## PLOT WITH DENDEXTEND
library(dendextend)

# Convert hclust to dendrogram
dend <- as.dendrogram(hc)

# Customize appearance
dend <- dend %>%
  set("branches_k_color", k = 3) %>%  # color branches by cluster
  set("branches_lwd", 3) %>%
  set("labels_cex", 1)

# Plot
plot(dend, main = "Fused Dendrogram of the Graphs")




## PLOT WITH GGDENDRO
library(ggdendro)

# Convert to dendrogram and then to ggplot
ggd <- as.dendrogram(hc)
ggd_data <- dendro_data(ggd)

# Plot
library(ggplot2)
ggplot(segment(ggd_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_minimal() +
  labs(title = "Dendrogram (ggplot2)", x = "", y = "Height")