

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

plots <- wrap_plots(
    p1, 
    p2, 
    # p3, 
    nrow = 1)


# === gene overlap structure between donors ===

euler1 <- per_donor_data[["logical_matrix_A549separate"]] %>%
  select(contains("A549")) %>%
  filter(!if_all(everything(), ~ .x == 0)) %>%
  plot_euler(
    returnAsFunction = T, 
    fillColors = rep("white", 10), showLabels = F, 
    plotTitle = "A549 replicates",
    aspectRatio = 0.5
    )

euler2 <- per_donor_data[["logical_matrix_A549separate"]] %>%
  select(contains("ALI")) %>%
  filter(!if_all(everything(), ~ .x == 0)) %>%
  plot_euler(
    returnAsFunction = T, 
    fillColors = rep("#D9D9D9", 10), 
    showLabels = F, 
    plotTitle = "ALI donors",
    aspectRatio = 0.5
    )

euler3 <- per_donor_data[["logical_matrix_A549separate"]] %>%
  select(contains("HBE")) %>%
  filter(!if_all(everything(), ~ .x == 0)) %>%
  plot_euler(
    returnAsFunction = T, 
    fillColors = rep("#A6CEE3", 10), 
    showLabels = F, 
    plotTitle = "HBE donors",
    aspectRatio = 0.5
    )

euler <- wrap_elements(euler1()) / wrap_elements(euler2()) / wrap_elements(euler3()) +
    plot_annotation(title = "C. overlap structure across celltypes")


# === heatmap of union of all 14621 genes regulated in any donor ===


colsplits = c(
  rep("A549", 4),
  rep("ALI", 5),
  rep("HBE", 5)
)

colors = circlize::colorRamp2(
  c(
    seq(
      quantile(clustered_matrices$all_donors_union, 0.02),
      -0.05,
      length = 75
    ),
    seq(
      -0.49,
      0.49,
      length = 50
    ),
    seq(
      0.5,
      quantile(clustered_matrices$all_donors_union, 0.98),
      length = 75
    )
  ),
  colorRampPalette(c("steelblue4", "steelblue2" , "white", "firebrick2", "firebrick4"))(200)
)


hm <- draw(
  Heatmap(
    clustered_matrices$all_donors_union,
    cluster_columns = F,
    cluster_rows = F,
    clustering_method_rows = "ward.D2",
    use_raster = T,
    show_row_names = F,
    column_split = colsplits,
    border = "black",
    col = colors,
    name = "log2fold"),
  column_title = "D. union of all IL1B-regulated genes"
)

temp_figs <- list(
  plots,
  euler,
  hm
)





