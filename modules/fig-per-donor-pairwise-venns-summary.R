df <- per_donor_data[["overlaps_summary"]] %>%
  filter(name != "A549 mean") %>%
  select(criterion, name, union, overlap, jaccard) %>%
  pivot_longer(cols = c(union, overlap, jaccard), names_to = "stat", values_to = "value") %>%
  mutate(
    group = if_else(criterion == "reg", "all regulated genes", "split by direction"),
    stat = factor(stat, levels = c("union", "overlap", "jaccard"))
  )

print(
  ggplot(df,
         aes(x = criterion, y = value)
  ) +
    facet_grid(stat~group, scales = "free") +
    stat_summary(geom = "bar", fun = "mean", width = 0.5, color = "black") +
    # stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.3) +
    geom_jitter(width = 0.05) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(y = NULL, x = NULL)
)

# t.test(
#   per_donor_data[["overlaps_summary"]] %>% 
#     filter(criterion == "upreg", name != "A549 mean") %>% 
#     pull(jaccard),
#     per_donor_data[["overlaps_summary"]] %>% 
#     filter(criterion == "downreg", name != "A549 mean") %>% 
#     pull(jaccard)
# )