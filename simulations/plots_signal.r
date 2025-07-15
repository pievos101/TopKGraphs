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
L_melt$signal = factor(L_melt$signal, 
		labels=c("0.01","0.05","0.10","0.20","0.30"))

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
L_melt$signal = factor(L_melt$signal, 
		labels=c("0.10","0.20","0.30","0.40","0.50"))

p1 = ggplot(L_melt, aes(x = signal, y = value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, outlier.fill = "white", outlier.color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Clustering performance",
    y = "Adjusted R-Index (ARI)",
    x = "Intra-cluster edge probability"
  ) #+
  #theme(legend.position = "none")

print(p1)