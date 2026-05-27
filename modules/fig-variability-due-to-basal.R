# === gene overlap structure between donors ====================================

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

print(
  wrap_elements(euler1()) /
    wrap_elements(euler2()) / 
    wrap_elements(euler3()) +
    plot_annotation(title = "A. overlap structure across cell type replicates")
)

# === heatmap of union of all 14621 genes regulated in any donor ===============

colsplits = c(
  rep("A549", 4),
  rep("ALI", 5),
  rep("HBE", 5)
)


colors = generic_l2f_heatmap_colors(
  clustered_matrices$all_donors_union
)

hm1 <- Heatmap(
  clustered_matrices$all_donors_union,
  cluster_columns = F,
  cluster_rows = F,
  use_raster = T,
  show_row_names = F,
  column_split = colsplits,
  border = "black",
  col = colors,
  name = "log2fold IL1B",
  width = unit(2.25, "inch")
)

colors = generic_l2f_heatmap_colors(
  clustered_matrices$all_donors_union_comparativeTPMs,
  colorscale = viridis::viridis(5)
)

hm2 <- Heatmap(
  clustered_matrices$all_donors_union_comparativeTPMs,
  cluster_columns = F,
  cluster_rows = F,
  use_raster = T,
  show_row_names = F,
  column_split = colsplits,
  border = "black",
  col = colors,
  name = "log2tpm at NS",
  width = unit(2.25, "inch")
)

draw(
  hm1 + hm2,
  column_title = "B. Genes regulated in 1+ samples"
)


# === order A549-regulated genes by standard deviation and check NS tpm ========

genes <- per_donor_data$thresholded_log2fold_A549separate %>%
  filter(
    celltype == "A549",
    meets_threshold
  ) %>%
  pull(Gene)

order <- per_donor_data$spread_A549separate %>%
  filter(
    Gene %in% genes,
    celltype == "A549"
  ) %>%
  arrange(sd) %>%
  pull(Gene)

l2f_mat <- per_donor_data$thresholded_log2fold_A549separate %>%
  filter(
    Gene %in% genes,
    celltype == "A549"
  ) %>%
  mutate(name = paste(celltype, rep)) %>%
  select(Gene, name, log2fold) %>%
  pivot_wider(names_from = "name", values_from = "log2fold") %>%
  column_to_rownames("Gene") %>%
  as.matrix() %>%
  .[order, ]

l2f_colors = generic_l2f_heatmap_colors(
  l2f_mat
)

tpm_mat <- tpm %>%
  filter(
    Gene %in% genes,
    treatment == "NS",
    celltype == "A549"
  ) %>%
  mutate(name = paste(celltype, rep)) %>%
  select(Gene, name, log2tpm) %>%
  pivot_wider(names_from = "name", values_from = "log2tpm") %>%
  column_to_rownames("Gene") %>%
  as.matrix() %>%
  .[order, ]

tpm_colors = generic_l2f_heatmap_colors(
  tpm_mat,
  viridis::viridis(5)
)

hm1 <- Heatmap(
  l2f_mat,
  cluster_columns = F,
  cluster_rows = F,
  use_raster = T,
  show_row_names = F,
  # column_split = colsplits,
  border = "black",
  col = l2f_colors,
  name = "log2fold IL1B",
  width = unit(0.75, "inch"),
  height = unit(4, "inch")
)

hm2 <- Heatmap(
  tpm_mat,
  cluster_columns = F,
  cluster_rows = F,
  use_raster = T,
  show_row_names = F,
  # column_split = colsplits,
  border = "black",
  col = tpm_colors,
  name = "log2tpm at NS",
  width = unit(0.75, "inch"),
  height = unit(4, "inch")
)

draw(
  hm1 + hm2,
  column_title = "C. A549 regulated genes ordered by SD"
)

# === sd vs basal tpm ==========================================================

genes <- per_donor_data$thresholded_log2fold_A549separate %>%
  filter(
    celltype == "A549",
    meets_threshold
  ) %>%
  pull(Gene)

tpm_df <- tpm %>%
  filter(
    celltype == "A549",
    treatment == "NS",
    Gene %in% genes
  ) %>%
  group_by(Gene, log2tpm) %>%
  summarize(log2tpm = mean(log2tpm), .groups = "drop") %>%
  filter(round(log2tpm, 3) != -3.322)

sd_df <- per_donor_data$spread_A549separate %>%
  filter(
    celltype == "A549",
    Gene %in% genes
  )

summary_df <- left_join(tpm_df, sd_df)

p1 <- ggplot(
  summary_df,
  aes(x = ntile(sd, 10), y = sd)
) +
  geom_boxplot(
    aes(group = ntile(sd, 10)),
    outliers = F
  ) +
  theme_minimal() + 
  theme(aspect.ratio = 0.5) +
  labs(
    x = "bin",
    y = expression("SD of " * log[2] * "fold across reps"),
    title = "D. binning A549 IL1B-induced genes by variability"
  ) +
  scale_x_continuous(breaks = seq(1, 10))

p2 <- ggplot(
  summary_df,
  aes(x = ntile(sd, 10), y = log2tpm)
) +
  geom_boxplot(
    aes(group = ntile(sd, 10)),
    outliers = F
  ) +
  theme_minimal() + 
  theme(aspect.ratio = 0.5) +
  labs(
    x = "bin",
    y = expression("basal expression (" * log[2] * "tpm at NS)"),
    title = "E. highly variable genes are lowly expressed"
  )

print(
  wrap_plots(p1 / p2)
)





