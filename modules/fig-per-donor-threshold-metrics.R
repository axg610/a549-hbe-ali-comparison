

# === plots for general per-donor metrics ===

p1 <- ggplot(
  per_donor_data[["counts_A549separate"]], 
  aes(y = name, x = n, fill = celltype, label = n)
) +
  geom_col(color = "black") +
  geom_text(hjust = 1.1) +
  theme_minimal() +
  labs(x = "genes with |log2fold| > 1", y = NULL, title = "A. regulated genes") +
  scale_fill_manual(values=  c("#FFFFFF", "#D9D9D9", "#A6CEE3")) +
  guides(fill = "none")

p2 <- ggplot(
  per_donor_data[["spread_A549separate"]],
  aes(y = celltype, x = sd, fill = celltype)
) +
  geom_jitter(size = 0.1, height = 0.3, alpha = 0.5) +
  geom_violin(trim = TRUE, alpha = 0.7, color = "black") +
  geom_boxplot(width = 0.2, outliers = F, alpha = 0.7, color = "black") +
  theme_minimal() +
  labs(y = NULL, x = "SD of log2fold across reps", title = "B. variability between replicates") +
  coord_cartesian(xlim = c(0, 4)) +
  scale_fill_manual(values=  c("#FFFFFF", "#D9D9D9", "#A6CEE3")) +
  guides(fill = "none")

# p3 <- ggplot(
#   per_donor_data[["spread_A549separate"]],
#   aes(y = celltype, x = percent_meets_threshold, fill = celltype)
# ) +
#   stat_summary(
#     geom = "errorbar", fun.data = "mean_sdl",
#     fun.args = list(mult = 1), 
#     width = 0.5
#   ) +
#   stat_summary(
#     geom = "col", fun = "mean",
#     color = "black"
#   ) +
#   theme_minimal() +
#   scale_fill_manual(values=  c("#FFFFFF", "#D9D9D9", "#A6CEE3")) +
#   labs(y = NULL, x = "% reps with |log2fold| > 1", title = "C. 'agreement'")

print(
  wrap_plots(
    p1, 
    p2, 
    # p3, 
    nrow = 1)
)

# === heatmap of union of all 14621 genes regulated in any donor ===

colsplits = c(
  rep("A549", 1),
  rep("ALI", 5),
  rep("HBE", 5)
)

colors = generic_l2f_heatmap_colors(clustered_matrices$all_donors_union_A5mean)

draw(
  Heatmap(
    clustered_matrices$all_donors_union_A5mean,
    cluster_columns = F,
    cluster_rows = F,
    use_raster = T,
    show_row_names = F,
    show_column_names = F,
    column_split = colsplits,
    border = "black",
    col = colors,
    name = "log2fold",
    width = unit(2.25, "inch")
  ),
  column_title = "C. all IL1B-regulated genes across cell type replicates"
)



