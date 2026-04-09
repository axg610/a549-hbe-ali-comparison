# plots

p1 <- ggplot(
  per_donor_data$counts, 
  aes(y = name, x = n, fill = celltype, label = n)
) +
  geom_col() +
  geom_text(hjust = 1.1) +
  theme_minimal() +
  labs(x = "genes with |log2fold| > 1", y = NULL, title = "A") +
  guides(fill = "none")

p2 <- ggplot(
  per_donor_data$spread,
  aes(y = celltype, x = sd, fill = celltype)
) +
  geom_jitter(size = 0.1, height = 0.3, alpha = 0.7) +
  geom_violin(trim = TRUE, alpha = 0.7, color = "black") +
  geom_boxplot(width = 0.2, outliers = F, alpha = 0.7, color = "black") +
  theme_minimal() +
  labs(y = NULL, x = "SD of log2fold across donors", title = "B") +
  coord_cartesian(xlim = c(0, 3)) +
  guides(fill = "none")

p3 <- ggplot(
  per_donor_data$spread,
  aes(y = celltype, x = percent_meets_threshold, fill = celltype)
) +
  stat_summary(
    geom = "errorbar", fun.data = "mean_sdl",
    fun.args = list(mult = 1), 
    width = 0.5
  ) +
  stat_summary(
    geom = "col", fun = "mean",
    color = "black"
  ) +
  theme_minimal() +
  labs(y = NULL, x = "% donors with |log2fold| > 1", title = "C")

print(
  wrap_plots(p1, p2, p3, nrow = 1)
)