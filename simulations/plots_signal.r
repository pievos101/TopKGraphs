one = read.table("001_inter_kNN.txt")
two = read.table("005_inter_kNN.txt")
three = read.table("010_inter_kNN.txt")
four = read.table("020_inter_kNN.txt")
five = read.table("030_inter_kNN.txt")

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
    title = "Classification performance",
    y = "Accuracy",
    x = "Inter-cluster edge probability"
  ) #+
  #theme(legend.position = "none")

print(p1)

######################################################


one = read.table("10_intra_kNN.txt")
two = read.table("20_intra_kNN.txt")
three = read.table("30_intra_kNN.txt")
four = read.table("40_intra_kNN.txt")
five = read.table("50_intra_kNN.txt")

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
    title = "Classification performance",
    y = "Accuracy",
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