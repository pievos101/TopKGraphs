# TopKGraphs
library(TopKLists)
library(igraph)
library(fastcluster) 
library(cluster) 

source("/home/bastian/GitHub/TopKGraphs/simulations/sim.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_SIL.R")
source("/home/bastian/GitHub/TopKGraphs/R/calc_BINARY.R")

# Varianz 
my_var <- 1
VVV <- c(0.1,0.5,1,5,10,15,20)
this_method = "ward.D"
max.k = 10
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
