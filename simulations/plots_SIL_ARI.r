


# -------------------------------
# Plot 1
# -------------------------------

RESULT_df = read.table("BreastCancer_SIL.txt")

custom_colors <- c(
  "TopKGraphs" = "#F8766D",  # reddish
  "Node2Vec"   = "#B79F00",  # yellow‑brownish / military‑greenish
  "Jaccard"    = "#00ba7c",  # greenish
  "Dice"       = "#00BFC4",  # cyan/teal
  "PageRank"   = "#C77CFF",  # lila / purple
  "Laplacian"  = "#FDE725"   # yellowish
)
# Ensure Method factor is correct
RESULT_df$Method <- factor(RESULT_df$Method,
                           levels = c("TopKGraphs", "Node2Vec",
                                      "Jaccard", "Dice", "PageRank", "Laplacian"))

# Make Metric a factor with Silhouette first
RESULT_df$Metric <- factor(RESULT_df$Metric,
                           levels = c("Silhouette","CalinskiHarabasz",
                           "DaviesBouldin")) 
# replace OtherMetric1, OtherMetric2 with the other metrics in your dataset

# Plot
ggplot(RESULT_df, aes(x = factor(k), y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(
    #title = "Cluster quality metrics vs ground truth",
    x = "k (kNN graph)",
    y = "Metric value"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = custom_colors)



# ======================================================
# 8. Plot 2
# ======================================================
RESULT_df = read.table("BreastCancer_ARI.txt")

# Ensure Method factor is correct
RESULT_df$Method <- factor(RESULT_df$Method,
                           levels = c("TopKGraphs", "Node2Vec",
                                    "Jaccard", "Dice", "PageRank", "Laplacian"))

ggplot(RESULT_df, aes(x = k, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal(base_size = 14) +
   scale_fill_manual(values = custom_colors)+
  labs(
    #title = "Clustering performance on BreastCancer dataset",
    x = "k (kNN graph)",
    y = "Score"
  )