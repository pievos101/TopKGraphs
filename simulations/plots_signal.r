one = read.table("05_wl.txt")
two = read.table("10_wl.txt")
three = read.table("20_wl.txt")
four = read.table("50_wl.txt")
five = read.table("100_wl.txt")

L = list()
L[[1]] = one
L[[2]] = two
L[[3]] = three
L[[4]] = four
L[[5]] = five
library(reshape)
library(ggplot2)

L_melt = melt(L)
colnames(L_melt) = c("Method", "value", "signal")


# Reorder factor (includes Laplace even if missing)
L_melt$Method <- factor(L_melt$Method, 
levels = c("TopKGraphs", "Node2Vec", 
"Jaccard", "Dice", "PageRank", "Laplacian"))

keep_methods <- c("TopKGraphs", "Node2Vec")

L_melt <- L_melt[L_melt$Method %in% keep_methods, ]


custom_colors <- c(
  "TopKGraphs" = "#F8766D",  # reddish
  "Node2Vec"   = "#B79F00"  # yellow‑brownish / military‑greenish
  #"Jaccard"    = "#00ba7c",  # greenish
  #"Dice"       = "#00BFC4",  # cyan/teal
  #"PageRank"   = "#C77CFF",  # lila / purple
  #"Laplacian"  = "#FDE725"   # yellowish
)


L_melt$signal = factor(
  L_melt$signal,
  labels = c("5", "10", "20", "50", "100")
)

# Mean
mean_df <- aggregate(
  value ~ Method + signal,
  data = L_melt,
  FUN = mean
)

# Standard error
se_df <- aggregate(
  value ~ Method + signal,
  data = L_melt,
  FUN = function(x) sd(x) / sqrt(length(x))
)

# Merge mean and SE
summary_df <- merge(mean_df, se_df,
                    by = c("Method", "signal"),
                    suffixes = c("_mean", "_se"))

p1 <- ggplot(summary_df,
             aes(x = signal, y = value_mean,
                 group = Method, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = value_mean - value_se,
        ymax = value_mean + value_se),
    width = 0.15,
    linewidth = 0.7
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = custom_colors) + 
  labs(
    #title = "Clustering performance",
    y = "Adjusted R-Index",
    x = "Walk Length"
  )

print(p1)
##################################


one = read.table("03_mu_kNN.txt")
two = read.table("05_mu_kNN.txt")
three = read.table("10_mu_kNN.txt")
four = read.table("20_mu_kNN.txt")
five = read.table("30_mu_kNN.txt")

L = list()
L[[1]] = one
L[[2]] = two
L[[3]] = three
L[[4]] = four
L[[5]] = five

library(reshape)
library(ggplot2)

L_melt = melt(L)
colnames(L_melt) = c("Method","value","signal")
L_melt$signal = factor(L_melt$signal, 
		labels=c("0.03","0.05","0.10","0.20","0.30"))

p1 = ggplot(L_melt, aes(x = signal, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Classification performance",
    y = "Accuracy",
    x = "Fraction of inter-community edges"
  ) #+
  #theme(legend.position = "none")

print(p1)

#####

one = read.table("03_mu.txt")
two = read.table("05_mu.txt")
three = read.table("10_mu.txt")
four = read.table("20_mu.txt")
five = read.table("30_mu.txt")

L = list()
L[[1]] = one
L[[2]] = two
L[[3]] = three
L[[4]] = four
L[[5]] = five

library(reshape)
library(ggplot2)

L_melt = melt(L)
colnames(L_melt) = c("Method","value","signal")


# Reorder factor (includes Laplace even if missing)
L_melt$Method <- factor(L_melt$Method, 
levels = c("TopKGraphs", "Node2Vec", 
"Jaccard", "Dice", "PageRank", "Laplacian"))

custom_colors <- c(
  "TopKGraphs" = "#F8766D",  # reddish
  "Node2Vec"   = "#B79F00",  # yellow‑brownish / military‑greenish
  "Jaccard"    = "#00ba7c",  # greenish
  "Dice"       = "#00BFC4",  # cyan/teal
  "PageRank"   = "#C77CFF",  # lila / purple
  "Laplacian"  = "#FDE725"   # yellowish
)

L_melt$signal = factor(L_melt$signal, 
		labels=c("0.03","0.05","0.10","0.20","0.30"))

p1 = ggplot(L_melt, aes(x = signal, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  #scale_fill_brewer(palette = "Set2") +
  labs(
    #title = "Clustering performance",
    y = "Adjusted R-Index (ARI)",
    x = "Fraction of inter-community edges"
  ) #+
  #theme(legend.position = "none")

print(p1)

#####

one = read.table("40_inter.txt")
two = read.table("50_inter.txt")
three = read.table("60_inter.txt")
four = read.table("70_inter.txt")

L = list()
L[[1]] = one
L[[2]] = two
L[[3]] = three
L[[4]] = four

library(reshape)
library(ggplot2)

L_melt = melt(L)
colnames(L_melt) = c("Method","value","signal")
L_melt$signal = factor(L_melt$signal, 
		labels=c("0.40","0.50","0.60","0.70"))

p1 = ggplot(L_melt, aes(x = signal, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Clustering performance",
    y = "Adjusted R-Index (ARI)",
    x = "Inter-cluster edge probability"
  ) #+
  #theme(legend.position = "none")

print(p1)


###
one = read.table("01_inter.txt")
two = read.table("05_inter.txt")
three = read.table("10_inter.txt")
four = read.table("20_inter.txt")
five = read.table("30_inter.txt")

L = list()
L[[1]] = one
L[[2]] = two
L[[3]] = three
L[[4]] = four
L[[5]] = five

library(reshape)
library(ggplot2)

L_melt = melt(L)
colnames(L_melt) = c("Method","value","signal")

# Reorder factor (includes Laplace even if missing)
L_melt$Method <- factor(L_melt$Method, 
levels = c("TopKGraphs", "Node2Vec", 
"Jaccard", "Dice", "PageRank", "Laplacian"))

custom_colors <- c(
  "TopKGraphs" = "#F8766D",  # reddish
  "Node2Vec"   = "#B79F00",  # yellow‑brownish / military‑greenish
  "Jaccard"    = "#00ba7c",  # greenish
  "Dice"       = "#00BFC4",  # cyan/teal
  "PageRank"   = "#C77CFF",  # lila / purple
  "Laplacian"  = "#FDE725"   # yellowish
)


L_melt$signal = factor(L_melt$signal, 
		labels=c("0.01","0.05","0.10","0.20","0.30"))

p1 = ggplot(L_melt, aes(x = signal, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  #scale_fill_brewer(palette = "Set2") +
  labs(
    #title = "Clustering performance",
    y = "Adjusted R-Index (ARI)",
    x = "Inter-cluster edge probability"
  ) #+
  #theme(legend.position = "none")

print(p1)

######################################################


one = read.table("10_intra.txt")
two = read.table("20_intra.txt")
three = read.table("30_intra.txt")
four = read.table("40_intra.txt")
five = read.table("50_intra.txt")

L = list()
L[[1]] = one
L[[2]] = two
L[[3]] = three
L[[4]] = four
L[[5]] = five

library(reshape)
library(ggplot2)

L_melt = melt(L)
colnames(L_melt) = c("Method","value","signal")

# Reorder factor (includes Laplace even if missing)
L_melt$Method <- factor(L_melt$Method, 
levels = c("TopKGraphs", "Node2Vec", 
"Jaccard", "Dice", "PageRank", "Laplacian"))

custom_colors <- c(
  "TopKGraphs" = "#F8766D",  # reddish
  "Node2Vec"   = "#B79F00",  # yellow‑brownish / military‑greenish
  "Jaccard"    = "#00ba7c",  # greenish
  "Dice"       = "#00BFC4",  # cyan/teal
  "PageRank"   = "#C77CFF",  # lila / purple
  "Laplacian"  = "#FDE725"   # yellowish
)

L_melt$signal = factor(L_melt$signal, 
		labels=c("0.10","0.20","0.30","0.40","0.50"))

p1 = ggplot(L_melt, aes(x = signal, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = custom_colors) +
  labs(
    #title = "Clustering performance",
    y = "Adjusted R-Index (ARI)",
    x = "Intra-cluster edge probability"
  ) #+
  #theme(legend.position = "none")

print(p1)

###########################################################


one = read.table("WD_2.txt")
two = read.table("WD_3.txt")
three = read.table("WD_5.txt")
four = read.table("WD_10.txt")
five = read.table("WD_20.txt")
six = read.table("WD_30.txt")
seven = read.table("WD_50.txt")


L = list()
L[[1]] = one
L[[2]] = two
L[[3]] = three
L[[4]] = four
L[[5]] = five
L[[6]] = six
L[[7]] = seven


library(reshape)
library(ggplot2)

L_melt = melt(L)
colnames(L_melt) = c("Method","value","signal")
L_melt$signal = factor(L_melt$signal, 
		labels=c("2","3","5","10","20","30","50"))

p1 = ggplot(L_melt, aes(x = signal, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Clustering performance",
    y = "Adjusted R-Index (ARI)",
    x = "Walk Length"
  ) +
  theme(legend.position = "bottom")

print(p1)


###########################################################

one = read.table("30_intra_2NN.txt")[,c(1,5)]
two = read.table("30_intra_3NN.txt")[,c(1,5)]
three = read.table("30_intra_5NN.txt")[,c(1,5)]
four = read.table("30_intra_10NN.txt")[,c(1,5)]
five = read.table("30_intra_20NN.txt")[,c(1,5)]

L = list()
L[[1]] = one
L[[2]] = two
L[[3]] = three
L[[4]] = four
L[[5]] = five

library(reshape)
library(ggplot2)

L_melt = melt(L)
colnames(L_melt) = c("Method","value","signal")
L_melt$signal = factor(L_melt$signal, 
		labels=c("2","3","5","10","20"))

p1 = ggplot(L_melt, aes(x = signal, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Classification performance",
    y = "Accuracy",
    x = "Walk Length"
  ) #+
  #theme(legend.position = "none")

print(p1)