
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

# plot

plot_euler(reg, plotTitle = "(A) IL1B-regulated genes at 6h") %>% print()
plot_euler(upreg, plotTitle = "(B) IL1B-upregulated genes at 6h") %>% print()
plot_euler(downreg, plotTitle = "(C) IL1B-downregulated genes at 6h") %>% print()