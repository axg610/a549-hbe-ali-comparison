
# knitr::include_graphics("data/il1b-venn.png")

# create logical matrices

reg <- dea %>%
  filter(
    treatment == "IL1B",
    FDR <= 0.05,
    abs(log2fold) >= 1
  ) %>%
  group_by(celltype) %>%
  summarize(Gene = list(unique(Gene)), .groups = "drop") %>%
  unnest(Gene) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "celltype", values_from = "value", values_fill = 0) %>%
  column_to_rownames("Gene")

upreg <- dea %>%
  filter(
    treatment == "IL1B",
    FDR <= 0.05,
    log2fold >= 1
  ) %>%
  group_by(celltype) %>%
  summarize(Gene = list(unique(Gene)), .groups = "drop") %>%
  unnest(Gene) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "celltype", values_from = "value", values_fill = 0) %>%
  column_to_rownames("Gene")

downreg <- dea %>%
  filter(
    treatment == "IL1B",
    FDR <= 0.05,
    log2fold <= -1
  ) %>%
  group_by(celltype) %>%
  summarize(Gene = list(unique(Gene)), .groups = "drop") %>%
  unnest(Gene) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "celltype", values_from = "value", values_fill = 0) %>%
  column_to_rownames("Gene")

reg_fold_only <- dea %>%
  filter(
    treatment == "IL1B",
    abs(log2fold) >= 1
  ) %>%
  group_by(celltype) %>%
  summarize(Gene = list(unique(Gene)), .groups = "drop") %>%
  unnest(Gene) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "celltype", values_from = "value", values_fill = 0) %>%
  column_to_rownames("Gene")

reg_fdr_only <- dea %>%
  filter(
    treatment == "IL1B",
    FDR <= 0.05
  ) %>%
  group_by(celltype) %>%
  summarize(Gene = list(unique(Gene)), .groups = "drop") %>%
  unnest(Gene) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "celltype", values_from = "value", values_fill = 0) %>%
  column_to_rownames("Gene")

# plot

plot_euler(reg, plotTitle = "A. regulated") %>% print()
plot_euler(upreg, plotTitle = "B. upregulated") %>% print()
plot_euler(downreg, plotTitle = "C. downregulated") %>% print()

plot_euler(reg_fold_only, plotTitle = "D. regulated, fold threshold only") %>% print()
plot_euler(reg_fdr_only, plotTitle = "E. regulated, FDR threshold only") %>% print()




