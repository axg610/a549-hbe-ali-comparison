# il1b-regulated genes cloned from old venns
# might need to simplify + scale up later

genes_of_interest <- dea %>%
  filter(
    treatment == "IL1B",
    FDR <= 0.05,
    abs(log2fold) >= 1
  ) %>%
  distinct(Gene, celltype) %>%
  mutate(
    celltype = recode(
      celltype,
      A549 = "A",
      ALI  = "L",
      HBE  = "H"
    ),
    value = 1
  ) %>%
  pivot_wider(
    names_from = celltype,
    values_from = value,
    values_fill = 0
  ) %>%
  mutate(
    group = paste0(
      ifelse(A == 1, "A", ""),
      ifelse(L == 1, "L", ""),
      ifelse(H == 1, "H", "")
    )
  ) %>%
  group_by(group) %>%
  summarize(genes = list(Gene), .groups = "drop") %>%
  deframe()

# make matrices for each venn lobe
# technically this should be in modules/compute_heavy_matrices.R
# but i don't think worth movings

matrices <- list()

for(grp in names(genes_of_interest)){
  
  mat <- tpm %>%
    filter(Gene %in% genes_of_interest[[grp]]) %>%
    filter(!(treatment %in% c("NS"))) %>%
    mutate(
      celltype = factor(celltype, levels = c("A549", "ALI", "HBE")),
      treatment = factor(treatment, levels = c("NS", "GC", "IL1B", "combo")),
      name = paste(celltype, treatment, rep)
    ) %>%
    arrange(celltype, treatment) %>%
    select(Gene, name, log2fold) %>%
    pivot_wider(
      names_from = name,
      values_from = log2fold
    ) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  mat <- mat[
    hclust(dist(mat), method = "ward.D2")$order,
  ]
  
  matrices[[grp]] <- mat
  
}

# calculate heatmap colors and column splits

colorscale <- generic_l2f_heatmap_colors(do.call(rbind, matrices))

columnsplits <- factor(
  c(
    rep("A549 GC", 4),
    rep("A549 IL1B", 4),
    rep("A549 combo", 4),
    rep("ALI GC", 5),
    rep("ALI IL1B", 5),
    rep("ALI combo", 5),
    rep("HBE GC", 5),
    rep("HBE IL1B", 5),
    rep("HBE combo", 5)
  ),
  levels = c(
    "A549 GC", "A549 IL1B", "A549 combo",
    "ALI GC", "ALI IL1B", "ALI combo",
    "HBE GC", "HBE IL1B", "HBE combo"
  )
)

# make heatmaps for each venn lobe

heatmaps <- list()

for(grp in names(matrices)){
  
  heatmaps[[grp]] <- Heatmap(
    
    matrices[[grp]],
    name = paste0("log2fold (", grp, ")"),
    col = colorscale,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    show_row_names = FALSE,
    show_column_names = FALSE,
    
    column_split = columnsplits,
    column_title_rot = 45,
    
    row_title = grp,
    
    border = "black"
  )
  
}

# write out heatmaps to deploy on RNAseq explorer

# saveRDS(heatmaps, "data/rds/il1b_reg_genes_heatmap.rds")
saveRDS(matrices, "data/rds/il1b_reg_venn_matrices.rds")

# draw heatmaps for the knit

print(
  heatmaps[["A"]]
)

print(
  heatmaps[["L"]] %v% heatmaps[["H"]] %v% heatmaps[["AL"]] %v%
    heatmaps[["AH"]] %v% heatmaps[["LH"]] %v% heatmaps[["ALH"]]
)