one = read.table("10.txt")
two = read.table("20.txt")
three = read.table("30.txt")
four = read.table("40.txt")
five = read.table("50.txt")

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