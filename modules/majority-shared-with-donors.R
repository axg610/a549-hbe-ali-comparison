datasets_to_process <- list(
  "regulated"     = per_donor_data[["overlaps"]],
  "upregulated"   = per_donor_data[["overlaps_upreg"]],
  "downregulated" = per_donor_data[["overlaps_downreg"]]
)

results <- data.frame()


for(nm in names(datasets_to_process)){
  
  df <- datasets_to_process[[nm]]
  
  out <- df %>%
    select(name, genes, n_genes) %>%
    rowwise() %>%
    mutate(
      other_genes  = list(unique(unlist(df$genes[df$name != name]))),
      shared_genes = list(intersect(genes, other_genes)),
      unique_genes = list(setdiff(genes, other_genes)),
      n_other      = length(other_genes),
      n_shared     = length(shared_genes),
      n_unique     = length(unique_genes)
    ) %>%
    ungroup() %>%
    select(
      name, 
      genes, other_genes, shared_genes, unique_genes, 
      contains("n_")
    ) %>%
    mutate(
      regulation = nm
    )
  
  results = rbind(results, out)
}

summary <- results %>%
  select(name, regulation, contains("n_")) %>%
  pivot_longer(
    cols = contains("n_"),
    names_to = "group",
    values_to = "n"
  ) %>%
  mutate(
    description = case_when(
      group == "n_genes"  ~ "donor regulatory pool",
      group == "n_unique" ~ "unique to donor",
      group == "n_shared" ~ "shared w/ 1+ other(s)",
      group == "n_other"  ~ "unique to others"
    ),
    group = factor(
      group,
      levels = c("n_unique", "n_shared", "n_genes", "n_other")
    ),
    regulation = factor(
      regulation,
      levels = c("regulated", "upregulated", "downregulated")
    )
  )

print(
  ggplot(
    summary %>% 
      filter(
        group %in% c("n_unique", "n_shared"),
        name == "A549 mean"
      ),
    aes(x = group, y = n, label = n, fill = regulation)
  ) +
    facet_wrap(
      ~regulation,
      labeller = labeller(
        regulation = c(
          "regulated" = "A. regulated",
          "upregulated" = "B. upregulated",
          "downregulated" = "C.downregulated"
        )
      )) +
    stat_summary(
      geom = "col", fun = "mean",
      color = "black",
      width = 0.75
    ) +
    geom_text(
      vjust = -0.5
    ) +
    scale_x_discrete(
      labels = c(
        "unique", 
        "shared")
    ) +
    labs(x = NULL, y = "number of genes", title = "A549 vs donors") +
    theme_minimal() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = c("#00BA38", "#F8766D", "#619CFF")) +
    guides(fill = "none")
)

print(
  ggplot(
    summary %>% 
      filter(group %in% c("n_unique", "n_shared")),
    aes(x = group, y = n, fill = regulation)
  ) +
    facet_wrap(
      ~regulation,
      labeller = labeller(
        regulation = c(
          "regulated" = "A. regulated",
          "upregulated" = "B. upregulated",
          "downregulated" = "C.downregulated"
        )
      )) +
    stat_summary(
      geom = "col", fun = "mean",
      color = "black",
      width = 0.75
    ) +
    geom_point() +
    geom_line(
      aes(group = name),
      linetype = "dashed"
    ) +
    scale_x_discrete(
      labels = c(
        "unique", 
        "shared")
    ) +
    labs(x = NULL, y = "number of genes", title = "all cells vs all others") +
    theme_minimal() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = c("#00BA38", "#F8766D", "#619CFF")) +
    guides(fill = "none")
)




