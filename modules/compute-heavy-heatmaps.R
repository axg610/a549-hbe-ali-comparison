
# === the union of all 14621 genes regulated in any donor, hclust by row ===

union_genes <- per_donor_data$logical_matrix_A549separate %>%
  rownames()

all_donors_union <- tpm %>%
  filter(
    Gene %in% union_genes,
    treatment == "IL1B"
  ) %>%
  mutate(name = paste(celltype, rep)) %>%
  select(Gene, name, log2fold) %>%
  pivot_wider(names_from = "name", values_from = "log2fold") %>%
  column_to_rownames("Gene") %>%
  as.matrix() %>%
  .[hclust(dist(.), method = "ward.D2")$order, ]



# === write objects ===

clustered_matrices <- list()

clustered_matrices[["all_donors_union"]] <- all_donors_union
